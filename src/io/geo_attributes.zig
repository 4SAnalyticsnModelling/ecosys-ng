const std = @import("std");
const utils = @import("utils");
const parser = @import("input_parser");
const max_path_len = 1024;
///Lat-lon ranges and tile specifications
pub const LatLonRangeAndTileSpecs = struct {
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
pub const SimInit = struct {
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
pub const GeoAttr = struct {
    allocator: std.mem.Allocator = undefined,
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
    pub const Record = struct {
        morton: u64,
        ix: u32,
        iy: u32,
        lat_ud: i32,
        lon_ud: i32,
        elev: f32,
        matc: f32,
    };
    const ParsedRecord = struct {
        record: Record,
        tile_ix: usize,
        tile_iy: usize,
    };

    pub fn loadAndWriteBin(self: *GeoAttr, lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs, file_name: []const u8) !void {
        var line = try parser.readNextDataLine(self.reader);
        var tokens = parser.Tokens{};
        try tokens.tokenizeLine(line, 1, "geographical attributes file path", file_name, self.err_log);
        var geo_buf: [max_path_len]u8 = undefined;
        @memcpy(geo_buf[0..tokens.items[0].len], tokens.items[0]);
        const geo_attr_filename: []const u8 = geo_buf[0..tokens.items[0].len];
        try self.writeHeader();

        var max_ix: usize = 0;
        var max_iy: usize = 0;
        var saw_record = false;
        {
            var geo_reader: utils.FileReader = utils.FileReader{};
            try geo_reader.open(self.err_log, geo_attr_filename);
            defer geo_reader.close();
            geo_reader.reader();
            while (true) {
                line = try parser.readNextDataLine(geo_reader.buf_reader);
                if (std.mem.eql(u8, line, "EndOfStream")) break;
                const parsed = try self.parseRecord(line, lat_lon_rng_n_tile_specs, geo_attr_filename);
                saw_record = true;
                max_ix = @max(max_ix, @as(usize, parsed.record.ix));
                max_iy = @max(max_iy, @as(usize, parsed.record.iy));
            }
        }

        if (saw_record) {
            self.nx = max_ix + 1;
            self.ny = max_iy + 1;
        } else {
            self.nx = 0;
            self.ny = 0;
        }

        const tiles_x: usize = if (saw_record) (@divFloor(max_ix, lat_lon_rng_n_tile_specs.ntx) + 1) else 0;
        const tiles_y: usize = if (saw_record) (@divFloor(max_iy, lat_lon_rng_n_tile_specs.nty) + 1) else 0;
        const tile_count = tiles_x * tiles_y;
        var tile_counts = try self.allocator.alloc(usize, tile_count);
        defer self.allocator.free(tile_counts);
        @memset(tile_counts, 0);

        {
            var geo_reader: utils.FileReader = utils.FileReader{};
            try geo_reader.open(self.err_log, geo_attr_filename);
            defer geo_reader.close();
            geo_reader.reader();
            while (true) {
                line = try parser.readNextDataLine(geo_reader.buf_reader);
                if (std.mem.eql(u8, line, "EndOfStream")) break;
                const parsed = try self.parseRecord(line, lat_lon_rng_n_tile_specs, geo_attr_filename);
                const tile_id = parsed.tile_iy * tiles_x + parsed.tile_ix;
                tile_counts[tile_id] += 1;
            }
        }

        var tile_id: usize = 0;
        while (tile_id < tile_count) : (tile_id += 1) {
            const count = tile_counts[tile_id];
            if (count == 0) continue;
            const tile_ix = tile_id % tiles_x;
            const tile_iy = tile_id / tiles_x;

            var records = try self.allocator.alloc(Record, count);
            defer self.allocator.free(records);
            var filled: usize = 0;

            {
                var geo_reader: utils.FileReader = utils.FileReader{};
                try geo_reader.open(self.err_log, geo_attr_filename);
                defer geo_reader.close();
                geo_reader.reader();
                while (true) {
                    line = try parser.readNextDataLine(geo_reader.buf_reader);
                    if (std.mem.eql(u8, line, "EndOfStream")) break;
                    const parsed = try self.parseRecord(line, lat_lon_rng_n_tile_specs, geo_attr_filename);
                    if (parsed.tile_ix != tile_ix or parsed.tile_iy != tile_iy) continue;
                    records[filled] = parsed.record;
                    filled += 1;
                }
            }

            if (filled != count) return error.RecordCountMismatch;
            std.sort.block(Record, records[0..filled], {}, GeoAttr.lessThanMorton);
            for (records[0..filled]) |record| {
                try self.writeRecord(record);
            }
        }

        return;
    }
    fn parseRecord(self: *GeoAttr, line: []const u8, lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs, geo_attr_filename: []const u8) !ParsedRecord {
        var tokens = parser.Tokens{};
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

        const ix = utils.toXY(self.lon_ud, lat_lon_rng_n_tile_specs.lon_min_ud, lat_lon_rng_n_tile_specs.dlon_ud);
        const iy = utils.toXY(self.lat_ud, lat_lon_rng_n_tile_specs.lat_min_ud, lat_lon_rng_n_tile_specs.dlat_ud);
        const tile_ix = @divFloor(ix, lat_lon_rng_n_tile_specs.ntx);
        const tile_iy = @divFloor(iy, lat_lon_rng_n_tile_specs.nty);

        const record = Record{
            .morton = @intCast(utils.morton2D(ix, iy)),
            .ix = @as(u32, @intCast(ix)),
            .iy = @as(u32, @intCast(iy)),
            .lat_ud = self.lat_ud,
            .lon_ud = self.lon_ud,
            .elev = self.elevation,
            .matc = self.matc,
        };

        return ParsedRecord{
            .record = record,
            .tile_ix = tile_ix,
            .tile_iy = tile_iy,
        };
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
    fn writeRecord(self: *GeoAttr, record: Record) !void {
        try self.writeIntLittle(u64, record.morton);
        try self.writeIntLittle(u32, record.ix);
        try self.writeIntLittle(u32, record.iy);
        try self.writeIntLittle(i32, record.lat_ud);
        try self.writeIntLittle(i32, record.lon_ud);
        try self.writeIntLittle(u32, @as(u32, @bitCast(record.elev)));
        try self.writeIntLittle(u32, @as(u32, @bitCast(record.matc)));
    }
    fn lessThanMorton(_: void, a: Record, b: Record) bool {
        return a.morton < b.morton;
    }
};

pub const GeoAttrBinData = struct {
    allocator: std.mem.Allocator = undefined,
    err_log: *std.Io.Writer = undefined,
    lat_ud: []i32 = undefined,
    lon_ud: []i32 = undefined,
    elev: []f32 = undefined,
    matc: []f32 = undefined,

    pub fn deinit(self: *GeoAttrBinData) void {
        _ = self;
    }

    pub fn readBinToArraysWithArena(self: *GeoAttrBinData, geo_attr: *const GeoAttr, bin_path: []const u8) !void {
        const n = geo_attr.nx * geo_attr.ny;

        self.lat_ud = try self.allocator.alloc(i32, n);
        self.lon_ud = try self.allocator.alloc(i32, n);
        self.elev = try self.allocator.alloc(f32, n);
        self.matc = try self.allocator.alloc(f32, n);

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
        var rec_buf: [32]u8 = undefined;
        while (true) {
            bin_reader.buf_reader.readSliceAll(&rec_buf) catch |err| {
                if (err == error.EndOfStream) break;
                return err;
            };
            const ix = std.mem.readInt(u32, rec_buf[8..12], .little);
            const iy = std.mem.readInt(u32, rec_buf[12..16], .little);
            const lat_ud_ = std.mem.readInt(i32, rec_buf[16..20], .little);
            const lon_ud_ = std.mem.readInt(i32, rec_buf[20..24], .little);
            const elev_bits = std.mem.readInt(u32, rec_buf[24..28], .little);
            const mat_bits = std.mem.readInt(u32, rec_buf[28..32], .little);
            const flat_id: usize = @as(usize, ix) + @as(usize, iy) * geo_attr.nx;
            if (flat_id >= n) return error.OutOfBounds;
            self.lat_ud[flat_id] = lat_ud_;
            self.lon_ud[flat_id] = lon_ud_;
            self.elev[flat_id] = @as(f32, @bitCast(elev_bits));
            self.matc[flat_id] = @as(f32, @bitCast(mat_bits));
        }
    }

    pub const TileHandler = *const fn (tile_ix: usize, tile_iy: usize, records: []const GeoAttr.Record) anyerror!void;

    pub const TileBatch = struct {
        tile_ix: usize,
        tile_iy: usize,
        records: []const GeoAttr.Record,
    };

    pub const TileReader = struct {
        allocator: std.mem.Allocator,
        err_log: *std.Io.Writer,
        lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs,
        bin_reader: utils.FileReader,
        records: std.ArrayList(GeoAttr.Record),
        has_pending: bool = false,
        pending_record: GeoAttr.Record = undefined,
        pending_tile_ix: usize = 0,
        pending_tile_iy: usize = 0,
        eof: bool = false,

        pub fn init(
            allocator: std.mem.Allocator,
            err_log: *std.Io.Writer,
            lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs,
            bin_path: []const u8,
        ) !TileReader {
            var bin_reader: utils.FileReader = utils.FileReader{};
            try bin_reader.open(err_log, bin_path);
            bin_reader.reader();
            var buf_header: [12]u8 = undefined;
            try bin_reader.buf_reader.readSliceAll(&buf_header);
            if (!std.mem.eql(u8, buf_header[0..8], "geo_attr")) return error.BadHeader;
            const version = std.mem.readInt(u16, buf_header[8..10], .little);
            const flags = std.mem.readInt(u16, buf_header[10..12], .little);
            if (version != 1 or flags != 1) return error.BadHeader;

            return TileReader{
                .allocator = allocator,
                .err_log = err_log,
                .lat_lon_rng_n_tile_specs = lat_lon_rng_n_tile_specs,
                .bin_reader = bin_reader,
                .records = try std.ArrayList(GeoAttr.Record).initCapacity(allocator, 0),
            };
        }

        pub fn deinit(self: *TileReader) void {
            self.records.deinit(self.allocator);
            self.bin_reader.close();
        }

        fn readNextRecord(self: *TileReader) !?GeoAttr.Record {
            var rec_buf: [32]u8 = undefined;
            self.bin_reader.buf_reader.readSliceAll(&rec_buf) catch |err| {
                if (err == error.EndOfStream) return null;
                return err;
            };

            return GeoAttr.Record{
                .morton = std.mem.readInt(u64, rec_buf[0..8], .little),
                .ix = std.mem.readInt(u32, rec_buf[8..12], .little),
                .iy = std.mem.readInt(u32, rec_buf[12..16], .little),
                .lat_ud = std.mem.readInt(i32, rec_buf[16..20], .little),
                .lon_ud = std.mem.readInt(i32, rec_buf[20..24], .little),
                .elev = @as(f32, @bitCast(std.mem.readInt(u32, rec_buf[24..28], .little))),
                .matc = @as(f32, @bitCast(std.mem.readInt(u32, rec_buf[28..32], .little))),
            };
        }

        pub fn next(self: *TileReader) !?TileBatch {
            if (self.eof and !self.has_pending) return null;

            self.records.clearRetainingCapacity();

            var current_tile_ix: usize = 0;
            var current_tile_iy: usize = 0;
            if (self.has_pending) {
                current_tile_ix = self.pending_tile_ix;
                current_tile_iy = self.pending_tile_iy;
                try self.records.append(self.allocator, self.pending_record);
                self.has_pending = false;
            } else {
                const rec_opt = try self.readNextRecord();
                if (rec_opt == null) {
                    self.eof = true;
                    return null;
                }
                const rec = rec_opt.?;
                current_tile_ix = @divFloor(@as(usize, rec.ix), self.lat_lon_rng_n_tile_specs.ntx);
                current_tile_iy = @divFloor(@as(usize, rec.iy), self.lat_lon_rng_n_tile_specs.nty);
                try self.records.append(self.allocator, rec);
            }

            while (true) {
                const rec_opt = try self.readNextRecord();
                if (rec_opt == null) {
                    self.eof = true;
                    break;
                }
                const rec = rec_opt.?;
                const tile_ix = @divFloor(@as(usize, rec.ix), self.lat_lon_rng_n_tile_specs.ntx);
                const tile_iy = @divFloor(@as(usize, rec.iy), self.lat_lon_rng_n_tile_specs.nty);
                if (tile_ix != current_tile_ix or tile_iy != current_tile_iy) {
                    self.has_pending = true;
                    self.pending_record = rec;
                    self.pending_tile_ix = tile_ix;
                    self.pending_tile_iy = tile_iy;
                    break;
                }
                try self.records.append(self.allocator, rec);
            }

            return TileBatch{
                .tile_ix = current_tile_ix,
                .tile_iy = current_tile_iy,
                .records = self.records.items,
            };
        }
    };

    pub fn readBinByTile(
        self: *GeoAttrBinData,
        lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs,
        bin_path: []const u8,
        handler: TileHandler,
    ) !void {
        var bin_reader: utils.FileReader = utils.FileReader{};
        try bin_reader.open(self.err_log, bin_path);
        defer bin_reader.close();
        bin_reader.reader();
        var buf_header: [12]u8 = undefined;
        try bin_reader.buf_reader.readSliceAll(&buf_header);

        if (!std.mem.eql(u8, buf_header[0..8], "geo_attr")) return error.BadHeader;
        const version = std.mem.readInt(u16, buf_header[8..10], .little);
        const flags = std.mem.readInt(u16, buf_header[10..12], .little);
        if (version != 1 or flags != 1) return error.BadHeader;

        var records = try std.ArrayList(GeoAttr.Record).initCapacity(self.allocator, 0);
        defer records.deinit(self.allocator);

        var rec_buf: [32]u8 = undefined;
        var has_tile = false;
        var current_tile_ix: usize = 0;
        var current_tile_iy: usize = 0;

        while (true) {
            bin_reader.buf_reader.readSliceAll(&rec_buf) catch |err| {
                if (err == error.EndOfStream) break;
                return err;
            };

            const morton = std.mem.readInt(u64, rec_buf[0..8], .little);
            const ix_u32 = std.mem.readInt(u32, rec_buf[8..12], .little);
            const iy_u32 = std.mem.readInt(u32, rec_buf[12..16], .little);
            const lat_ud_ = std.mem.readInt(i32, rec_buf[16..20], .little);
            const lon_ud_ = std.mem.readInt(i32, rec_buf[20..24], .little);
            const elev_bits = std.mem.readInt(u32, rec_buf[24..28], .little);
            const mat_bits = std.mem.readInt(u32, rec_buf[28..32], .little);

            const ix = @as(usize, ix_u32);
            const iy = @as(usize, iy_u32);
            const tile_ix = @divFloor(ix, lat_lon_rng_n_tile_specs.ntx);
            const tile_iy = @divFloor(iy, lat_lon_rng_n_tile_specs.nty);

            if (!has_tile) {
                has_tile = true;
                current_tile_ix = tile_ix;
                current_tile_iy = tile_iy;
            } else if (tile_ix != current_tile_ix or tile_iy != current_tile_iy) {
                try handler(current_tile_ix, current_tile_iy, records.items);
                records.clearRetainingCapacity();
                current_tile_ix = tile_ix;
                current_tile_iy = tile_iy;
            }

            try records.append(self.allocator, GeoAttr.Record{
                .morton = morton,
                .ix = ix_u32,
                .iy = iy_u32,
                .lat_ud = lat_ud_,
                .lon_ud = lon_ud_,
                .elev = @as(f32, @bitCast(elev_bits)),
                .matc = @as(f32, @bitCast(mat_bits)),
            });
        }

        if (has_tile) {
            try handler(current_tile_ix, current_tile_iy, records.items);
        }
    }

    pub fn readBinByTileAndFillArrays(
        self: *GeoAttrBinData,
        geo_attr: *const GeoAttr,
        lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs,
        bin_path: []const u8,
        handler: TileHandler,
    ) !void {
        const n = geo_attr.nx * geo_attr.ny;
        self.lat_ud = try self.allocator.alloc(i32, n);
        self.lon_ud = try self.allocator.alloc(i32, n);
        self.elev = try self.allocator.alloc(f32, n);
        self.matc = try self.allocator.alloc(f32, n);

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

        if (!std.mem.eql(u8, buf_header[0..8], "geo_attr")) return error.BadHeader;
        const version = std.mem.readInt(u16, buf_header[8..10], .little);
        const flags = std.mem.readInt(u16, buf_header[10..12], .little);
        if (version != 1 or flags != 1) return error.BadHeader;

        var records = try std.ArrayList(GeoAttr.Record).initCapacity(self.allocator, 0);
        defer records.deinit(self.allocator);

        var rec_buf: [32]u8 = undefined;
        var has_tile = false;
        var current_tile_ix: usize = 0;
        var current_tile_iy: usize = 0;

        while (true) {
            bin_reader.buf_reader.readSliceAll(&rec_buf) catch |err| {
                if (err == error.EndOfStream) break;
                return err;
            };

            const morton = std.mem.readInt(u64, rec_buf[0..8], .little);
            const ix_u32 = std.mem.readInt(u32, rec_buf[8..12], .little);
            const iy_u32 = std.mem.readInt(u32, rec_buf[12..16], .little);
            const lat_ud_ = std.mem.readInt(i32, rec_buf[16..20], .little);
            const lon_ud_ = std.mem.readInt(i32, rec_buf[20..24], .little);
            const elev_bits = std.mem.readInt(u32, rec_buf[24..28], .little);
            const mat_bits = std.mem.readInt(u32, rec_buf[28..32], .little);

            const ix = @as(usize, ix_u32);
            const iy = @as(usize, iy_u32);
            const tile_ix = @divFloor(ix, lat_lon_rng_n_tile_specs.ntx);
            const tile_iy = @divFloor(iy, lat_lon_rng_n_tile_specs.nty);

            const flat_id: usize = ix + iy * geo_attr.nx;
            if (flat_id >= n) return error.OutOfBounds;
            self.lat_ud[flat_id] = lat_ud_;
            self.lon_ud[flat_id] = lon_ud_;
            self.elev[flat_id] = @as(f32, @bitCast(elev_bits));
            self.matc[flat_id] = @as(f32, @bitCast(mat_bits));

            if (!has_tile) {
                has_tile = true;
                current_tile_ix = tile_ix;
                current_tile_iy = tile_iy;
            } else if (tile_ix != current_tile_ix or tile_iy != current_tile_iy) {
                try handler(current_tile_ix, current_tile_iy, records.items);
                records.clearRetainingCapacity();
                current_tile_ix = tile_ix;
                current_tile_iy = tile_iy;
            }

            try records.append(self.allocator, GeoAttr.Record{
                .morton = morton,
                .ix = ix_u32,
                .iy = iy_u32,
                .lat_ud = lat_ud_,
                .lon_ud = lon_ud_,
                .elev = @as(f32, @bitCast(elev_bits)),
                .matc = @as(f32, @bitCast(mat_bits)),
            });
        }

        if (has_tile) {
            try handler(current_tile_ix, current_tile_iy, records.items);
        }
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
