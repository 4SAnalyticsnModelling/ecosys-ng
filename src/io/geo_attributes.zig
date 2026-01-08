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
        try parser.boundsCheck(error.TileSizeTooBig, .{self.ntx > 65536 or self.nty > 65536}, "tile specifications", file_name, self.err_log);
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
    nx: usize = 0, //max domain ix
    ny: usize = 0, //max domain iy
    pub const Record = struct {
        morton: u32,
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

        try parser.boundsCheck(error.TooManyGrids, .{self.nx > 65536 or self.ny > 65536}, "grid locations", geo_attr_filename, self.err_log);

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

            try parser.boundsCheck(error.RecordCountMismatch, .{filled != count}, "counting geospatial attributes records", geo_attr_filename, self.err_log);

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
            .morton = @as(u32, @intCast(utils.morton2D(ix, iy))),
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
    fn writeHeader(self: *GeoAttr) !void {
        try self.bin_writer.writeAll("geo_attr");
        try self.writeIntLittle(u16, 1); //version
        try self.writeIntLittle(u16, 1); //flags: little endian
    }
    fn writeIntLittle(self: *GeoAttr, comptime T: type, value: T) !void {
        var buf: [@sizeOf(T)]u8 = undefined;
        std.mem.writeInt(T, &buf, value, .little);
        try self.bin_writer.writeAll(&buf);
    }
    fn writeRecord(self: *GeoAttr, record: Record) !void {
        try self.writeIntLittle(u32, record.morton);
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

    pub const TileArrayHandler = *const fn (
        tile_ix: usize,
        tile_iy: usize,
        tile_nx: usize,
        tile_ny: usize,
        halo_stride: usize,
        tile_lat_ud: []const i32,
        tile_lon_ud: []const i32,
        tile_elev: []const f32,
        tile_matc: []const f32,
        records: []const GeoAttr.Record
    ) anyerror!void;

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
            bin_path: []const u8
        ) !TileReader {
            var bin_reader: utils.FileReader = utils.FileReader{};
            try bin_reader.open(err_log, bin_path);
            bin_reader.reader();
            var buf_header: [12]u8 = undefined;
            try bin_reader.buf_reader.readSliceAll(&buf_header);
            try parser.boundsCheck(error.BadHeader, .{!std.mem.eql(u8, buf_header[0..8], "geo_attr")}, "reading binary data header for geospatial attributes", bin_path, err_log);
            const version = std.mem.readInt(u16, buf_header[8..10], .little);
            const flags = std.mem.readInt(u16, buf_header[10..12], .little);
            try parser.boundsCheck(error.BadHeader, .{version != 1 or flags != 1}, "reading binary data header for geospatial attributes", bin_path, err_log);

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
            var rec_buf: [28]u8 = undefined;
            self.bin_reader.buf_reader.readSliceAll(&rec_buf) catch |err| {
                if (err == error.EndOfStream) return null;
                return err;
            };

            return GeoAttr.Record{
                .morton = std.mem.readInt(u32, rec_buf[0..4], .little),
                .ix = std.mem.readInt(u32, rec_buf[4..8], .little),
                .iy = std.mem.readInt(u32, rec_buf[8..12], .little),
                .lat_ud = std.mem.readInt(i32, rec_buf[12..16], .little),
                .lon_ud = std.mem.readInt(i32, rec_buf[16..20], .little),
                .elev = @as(f32, @bitCast(std.mem.readInt(u32, rec_buf[20..24], .little))),
                .matc = @as(f32, @bitCast(std.mem.readInt(u32, rec_buf[24..28], .little))),
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

    const TileRange = struct {
        start: usize,
        len: usize,
    };

    const RowBuffer = struct {
        allocator: std.mem.Allocator,
        tiles_x: usize,
        tile_stride: usize,
        tile_height: usize,
        tile_size: usize,
        lat_ud: []i32,
        lon_ud: []i32,
        elev: []f32,
        matc: []f32,
        tile_nx: []usize,
        tile_ny: []usize,
        present: []bool,
        tile_records: []TileRange,
        records: std.ArrayList(GeoAttr.Record),

        fn init(allocator: std.mem.Allocator, tiles_x: usize, tile_stride: usize, tile_height: usize) !RowBuffer {
            var row = RowBuffer{
                .allocator = allocator,
                .tiles_x = tiles_x,
                .tile_stride = tile_stride,
                .tile_height = tile_height,
                .tile_size = tile_stride * tile_height,
                .lat_ud = undefined,
                .lon_ud = undefined,
                .elev = undefined,
                .matc = undefined,
                .tile_nx = undefined,
                .tile_ny = undefined,
                .present = undefined,
                .tile_records = undefined,
                .records = undefined,
            };
            const row_size = tiles_x * row.tile_size;
            row.lat_ud = try allocator.alloc(i32, row_size);
            errdefer allocator.free(row.lat_ud);
            row.lon_ud = try allocator.alloc(i32, row_size);
            errdefer allocator.free(row.lon_ud);
            row.elev = try allocator.alloc(f32, row_size);
            errdefer allocator.free(row.elev);
            row.matc = try allocator.alloc(f32, row_size);
            errdefer allocator.free(row.matc);
            row.tile_nx = try allocator.alloc(usize, tiles_x);
            errdefer allocator.free(row.tile_nx);
            row.tile_ny = try allocator.alloc(usize, tiles_x);
            errdefer allocator.free(row.tile_ny);
            row.present = try allocator.alloc(bool, tiles_x);
            errdefer allocator.free(row.present);
            row.tile_records = try allocator.alloc(TileRange, tiles_x);
            errdefer allocator.free(row.tile_records);
            row.records = try std.ArrayList(GeoAttr.Record).initCapacity(allocator, 0);
            return row;
        }

        fn deinit(self: *RowBuffer) void {
            self.records.deinit(self.allocator);
            self.allocator.free(self.tile_records);
            self.allocator.free(self.present);
            self.allocator.free(self.tile_ny);
            self.allocator.free(self.tile_nx);
            self.allocator.free(self.matc);
            self.allocator.free(self.elev);
            self.allocator.free(self.lon_ud);
            self.allocator.free(self.lat_ud);
        }

        fn reset(self: *RowBuffer) void {
            @memset(self.present, false);
            self.records.clearRetainingCapacity();
        }

        fn tileBase(self: *const RowBuffer, tile_ix: usize) usize {
            return tile_ix * self.tile_size;
        }

        fn clearTileArrays(
            lat_buf: []i32,
            lon_buf: []i32,
            elev_buf: []f32,
            mat_buf: []f32,
            lat_sentinel: i32,
            lon_sentinel: i32,
            elev_nan: f32,
            mat_nan: f32
        ) void {
            for (lat_buf, lon_buf, elev_buf, mat_buf) |*lat_p, *lon_p, *elev_p, *mat_p| {
                lat_p.* = lat_sentinel;
                lon_p.* = lon_sentinel;
                elev_p.* = elev_nan;
                mat_p.* = mat_nan;
            }
        }

        fn beginTile(
            self: *RowBuffer,
            tile_ix: usize,
            tile_nx: usize,
            tile_ny: usize,
            lat_sentinel: i32,
            lon_sentinel: i32,
            elev_nan: f32,
            mat_nan: f32
        ) void {
            self.present[tile_ix] = true;
            self.tile_nx[tile_ix] = tile_nx;
            self.tile_ny[tile_ix] = tile_ny;
            self.tile_records[tile_ix] = TileRange{
                .start = self.records.items.len,
                .len = 0,
            };
            const base = self.tileBase(tile_ix);
            const tile_end = base + self.tile_size;
            RowBuffer.clearTileArrays(
                self.lat_ud[base..tile_end],
                self.lon_ud[base..tile_end],
                self.elev[base..tile_end],
                self.matc[base..tile_end],
                lat_sentinel,
                lon_sentinel,
                elev_nan,
                mat_nan
            );
        }

        fn endTile(self: *RowBuffer, tile_ix: usize) void {
            const range = &self.tile_records[tile_ix];
            range.len = self.records.items.len - range.start;
        }
    };

    pub fn readBinByTileAndFillArrays(
        self: *GeoAttrBinData,
        geo_attr: *const GeoAttr,
        lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs,
        bin_path: []const u8,
        handler: TileArrayHandler
    ) !void {
        const tile_stride = lat_lon_rng_n_tile_specs.ntx;
        const tile_height = lat_lon_rng_n_tile_specs.nty;
        const tile_size = tile_stride * tile_height;
        const halo_stride = tile_stride + 2;
        const halo_height = tile_height + 2;
        const halo_size = halo_stride * halo_height;

        self.lat_ud = try self.allocator.alloc(i32, halo_size);
        self.lon_ud = try self.allocator.alloc(i32, halo_size);
        self.elev = try self.allocator.alloc(f32, halo_size);
        self.matc = try self.allocator.alloc(f32, halo_size);

        const lat_sentinel = std.math.minInt(i32);
        const lon_sentinel = std.math.minInt(i32);
        const elev_nan = std.math.nan(f32);
        const mat_nan = std.math.nan(f32);

        const tiles_x: usize = if (geo_attr.nx > 0 and tile_stride > 0) (@divFloor(geo_attr.nx - 1, tile_stride) + 1) else 0;

        var row_prev = try RowBuffer.init(self.allocator, tiles_x, tile_stride, tile_height);
        defer row_prev.deinit();
        var row_curr = try RowBuffer.init(self.allocator, tiles_x, tile_stride, tile_height);
        defer row_curr.deinit();
        var row_next = try RowBuffer.init(self.allocator, tiles_x, tile_stride, tile_height);
        defer row_next.deinit();

        row_prev.reset();
        row_curr.reset();
        row_next.reset();

        const Helpers = struct {
            fn clearHalo(
                lat_buf: []i32,
                lon_buf: []i32,
                elev_buf: []f32,
                mat_buf: []f32,
                halo_lat_sentinel: i32,
                halo_lon_sentinel: i32,
                halo_elev_nan: f32,
                halo_mat_nan: f32
            ) void {
                for (lat_buf, lon_buf, elev_buf, mat_buf) |*lat_p, *lon_p, *elev_p, *mat_p| {
                    lat_p.* = halo_lat_sentinel;
                    lon_p.* = halo_lon_sentinel;
                    elev_p.* = halo_elev_nan;
                    mat_p.* = halo_mat_nan;
                }
            }

            fn processRow(
                handler_cb: TileArrayHandler,
                tile_stride_local: usize,
                halo_stride_local: usize,
                tiles_x_count: usize,
                halo_lat_sentinel: i32,
                halo_lon_sentinel: i32,
                halo_elev_nan: f32,
                halo_mat_nan: f32,
                row_prev_buf: ?*RowBuffer,
                row_curr_buf: *RowBuffer,
                row_next_buf: ?*RowBuffer,
                row_curr_iy: usize,
                halo_lat: []i32,
                halo_lon: []i32,
                halo_elev: []f32,
                halo_mat: []f32
            ) !void {
                var tile_ix: usize = 0;
                while (tile_ix < tiles_x_count) : (tile_ix += 1) {
                    if (!row_curr_buf.present[tile_ix]) continue;

                    const tile_nx = row_curr_buf.tile_nx[tile_ix];
                    const tile_ny = row_curr_buf.tile_ny[tile_ix];
                    if (tile_nx == 0 or tile_ny == 0) continue;

                    clearHalo(halo_lat, halo_lon, halo_elev, halo_mat, halo_lat_sentinel, halo_lon_sentinel, halo_elev_nan, halo_mat_nan);

                    const base = row_curr_buf.tileBase(tile_ix);
                    var y: usize = 0;
                    while (y < tile_ny) : (y += 1) {
                        var x: usize = 0;
                        while (x < tile_nx) : (x += 1) {
                            const src_idx = base + x + y * tile_stride_local;
                            const dst_idx = (x + 1) + (y + 1) * halo_stride_local;
                            halo_lat[dst_idx] = row_curr_buf.lat_ud[src_idx];
                            halo_lon[dst_idx] = row_curr_buf.lon_ud[src_idx];
                            halo_elev[dst_idx] = row_curr_buf.elev[src_idx];
                            halo_mat[dst_idx] = row_curr_buf.matc[src_idx];
                        }
                    }

                    if (tile_ix > 0 and row_curr_buf.present[tile_ix - 1]) {
                        const neigh_ix = tile_ix - 1;
                        const neigh_nx = row_curr_buf.tile_nx[neigh_ix];
                        const neigh_ny = row_curr_buf.tile_ny[neigh_ix];
                        const copy_ny = @min(tile_ny, neigh_ny);
                        const neigh_base = row_curr_buf.tileBase(neigh_ix);
                        var wy: usize = 0;
                        while (wy < copy_ny) : (wy += 1) {
                            const src_idx = neigh_base + (neigh_nx - 1) + wy * tile_stride_local;
                            const dst_idx = 0 + (wy + 1) * halo_stride_local;
                            halo_lat[dst_idx] = row_curr_buf.lat_ud[src_idx];
                            halo_lon[dst_idx] = row_curr_buf.lon_ud[src_idx];
                            halo_elev[dst_idx] = row_curr_buf.elev[src_idx];
                            halo_mat[dst_idx] = row_curr_buf.matc[src_idx];
                        }
                    }

                    if (tile_ix + 1 < tiles_x_count and row_curr_buf.present[tile_ix + 1]) {
                        const neigh_ix = tile_ix + 1;
                        const neigh_ny = row_curr_buf.tile_ny[neigh_ix];
                        const copy_ny = @min(tile_ny, neigh_ny);
                        const neigh_base = row_curr_buf.tileBase(neigh_ix);
                        var ey: usize = 0;
                        while (ey < copy_ny) : (ey += 1) {
                            const src_idx = neigh_base + 0 + ey * tile_stride_local;
                            const dst_idx = (tile_nx + 1) + (ey + 1) * halo_stride_local;
                            halo_lat[dst_idx] = row_curr_buf.lat_ud[src_idx];
                            halo_lon[dst_idx] = row_curr_buf.lon_ud[src_idx];
                            halo_elev[dst_idx] = row_curr_buf.elev[src_idx];
                            halo_mat[dst_idx] = row_curr_buf.matc[src_idx];
                        }
                    }

                    if (row_prev_buf) |prev| {
                        if (prev.present[tile_ix]) {
                            const neigh_nx = prev.tile_nx[tile_ix];
                            const neigh_ny = prev.tile_ny[tile_ix];
                            const copy_nx = @min(tile_nx, neigh_nx);
                            const neigh_base = prev.tileBase(tile_ix);
                            const row_y = neigh_ny - 1;
                            var sx: usize = 0;
                            while (sx < copy_nx) : (sx += 1) {
                                const src_idx = neigh_base + sx + row_y * tile_stride_local;
                                const dst_idx = (sx + 1) + 0 * halo_stride_local;
                                halo_lat[dst_idx] = prev.lat_ud[src_idx];
                                halo_lon[dst_idx] = prev.lon_ud[src_idx];
                                halo_elev[dst_idx] = prev.elev[src_idx];
                                halo_mat[dst_idx] = prev.matc[src_idx];
                            }
                        }
                    }

                    if (row_next_buf) |next| {
                        if (next.present[tile_ix]) {
                            const neigh_nx = next.tile_nx[tile_ix];
                            const copy_nx = @min(tile_nx, neigh_nx);
                            const neigh_base = next.tileBase(tile_ix);
                            var nx: usize = 0;
                            while (nx < copy_nx) : (nx += 1) {
                                const src_idx = neigh_base + nx + 0 * tile_stride_local;
                                const dst_idx = (nx + 1) + (tile_ny + 1) * halo_stride_local;
                                halo_lat[dst_idx] = next.lat_ud[src_idx];
                                halo_lon[dst_idx] = next.lon_ud[src_idx];
                                halo_elev[dst_idx] = next.elev[src_idx];
                                halo_mat[dst_idx] = next.matc[src_idx];
                            }
                        }
                    }

                    const range = row_curr_buf.tile_records[tile_ix];
                    const records_slice = row_curr_buf.records.items[range.start .. range.start + range.len];
                    try handler_cb(
                        tile_ix,
                        row_curr_iy,
                        tile_nx,
                        tile_ny,
                        halo_stride_local,
                        halo_lat,
                        halo_lon,
                        halo_elev,
                        halo_mat,
                        records_slice
                    );
                }
            }
        };

        var bin_reader: utils.FileReader = utils.FileReader{};
        try bin_reader.open(self.err_log, bin_path);
        defer bin_reader.close();
        bin_reader.reader();
        var buf_header: [12]u8 = undefined;
        try bin_reader.buf_reader.readSliceAll(&buf_header);

        try parser.boundsCheck(error.BadHeader, .{!std.mem.eql(u8, buf_header[0..8], "geo_attr")}, "reading binary data header for geospatial attributes", bin_path, self.err_log);
        const version = std.mem.readInt(u16, buf_header[8..10], .little);
        const flags = std.mem.readInt(u16, buf_header[10..12], .little);
        try parser.boundsCheck(error.BadHeader, .{version != 1 or flags != 1}, "reading binary data header for geospatial attributes", bin_path, self.err_log);

        var rec_buf: [28]u8 = undefined;
        var has_tile = false;
        var current_tile_ix: usize = 0;
        var current_tile_iy: usize = 0;
        var row_prev_iy: ?usize = null;
        var row_curr_iy: ?usize = null;
        var row_next_iy: ?usize = null;

        while (true) {
            bin_reader.buf_reader.readSliceAll(&rec_buf) catch |err| {
                if (err == error.EndOfStream) break;
                return err;
            };

            const morton = std.mem.readInt(u32, rec_buf[0..4], .little);
            const ix_u32 = std.mem.readInt(u32, rec_buf[4..8], .little);
            const iy_u32 = std.mem.readInt(u32, rec_buf[8..12], .little);
            const lat_ud_ = std.mem.readInt(i32, rec_buf[12..16], .little);
            const lon_ud_ = std.mem.readInt(i32, rec_buf[16..20], .little);
            const elev_bits = std.mem.readInt(u32, rec_buf[20..24], .little);
            const mat_bits = std.mem.readInt(u32, rec_buf[24..28], .little);

            const ix = @as(usize, ix_u32);
            const iy = @as(usize, iy_u32);
            const tile_ix = @divFloor(ix, tile_stride);
            const tile_iy = @divFloor(iy, tile_height);

            if (!has_tile) {
                has_tile = true;
                current_tile_ix = tile_ix;
                current_tile_iy = tile_iy;
                if (row_next_iy == null) {
                    row_next_iy = tile_iy;
                    row_next.reset();
                }
                const base_x = tile_ix * tile_stride;
                const base_y = tile_iy * tile_height;
                const remaining_x = if (geo_attr.nx > base_x) geo_attr.nx - base_x else 0;
                const remaining_y = if (geo_attr.ny > base_y) geo_attr.ny - base_y else 0;
                const tile_nx = @min(tile_stride, remaining_x);
                const tile_ny = @min(tile_height, remaining_y);
                row_next.beginTile(tile_ix, tile_nx, tile_ny, lat_sentinel, lon_sentinel, elev_nan, mat_nan);
            } else if (tile_ix != current_tile_ix or tile_iy != current_tile_iy) {
                row_next.endTile(current_tile_ix);

                if (tile_iy != current_tile_iy) {
                    if (row_curr_iy == null) {
                        std.mem.swap(RowBuffer, &row_curr, &row_next);
                        row_curr_iy = row_next_iy;
                        row_next_iy = null;
                        row_next.reset();
                    } else {
                        try Helpers.processRow(
                            handler,
                            tile_stride,
                            halo_stride,
                            tiles_x,
                            lat_sentinel,
                            lon_sentinel,
                            elev_nan,
                            mat_nan,
                            if (row_prev_iy != null) &row_prev else null,
                            &row_curr,
                            &row_next,
                            row_curr_iy.?,
                            self.lat_ud,
                            self.lon_ud,
                            self.elev,
                            self.matc
                        );

                        const tmp = row_prev;
                        row_prev = row_curr;
                        row_curr = row_next;
                        row_next = tmp;
                        row_prev_iy = row_curr_iy;
                        row_curr_iy = row_next_iy;
                        row_next_iy = null;
                        row_next.reset();
                    }
                }

                if (row_next_iy == null) {
                    row_next_iy = tile_iy;
                    row_next.reset();
                }
                const base_x = tile_ix * tile_stride;
                const base_y = tile_iy * tile_height;
                const remaining_x = if (geo_attr.nx > base_x) geo_attr.nx - base_x else 0;
                const remaining_y = if (geo_attr.ny > base_y) geo_attr.ny - base_y else 0;
                const tile_nx = @min(tile_stride, remaining_x);
                const tile_ny = @min(tile_height, remaining_y);
                row_next.beginTile(tile_ix, tile_nx, tile_ny, lat_sentinel, lon_sentinel, elev_nan, mat_nan);
                current_tile_ix = tile_ix;
                current_tile_iy = tile_iy;
            }

            const local_ix = ix - tile_ix * tile_stride;
            const local_iy = iy - tile_iy * tile_height;
            const flat_id: usize = local_ix + local_iy * tile_stride;
            if (flat_id >= tile_size) return error.OutOfBounds;
            const base = row_next.tileBase(tile_ix);
            const tile_idx = base + flat_id;
            row_next.lat_ud[tile_idx] = lat_ud_;
            row_next.lon_ud[tile_idx] = lon_ud_;
            row_next.elev[tile_idx] = @as(f32, @bitCast(elev_bits));
            row_next.matc[tile_idx] = @as(f32, @bitCast(mat_bits));

            try row_next.records.append(self.allocator, GeoAttr.Record{
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
            row_next.endTile(current_tile_ix);
        }

        if (row_next_iy != null) {
            if (row_curr_iy == null) {
                std.mem.swap(RowBuffer, &row_curr, &row_next);
                row_curr_iy = row_next_iy;
                row_next_iy = null;
                row_next.reset();
            } else {
                try Helpers.processRow(
                    handler,
                    tile_stride,
                    halo_stride,
                    tiles_x,
                    lat_sentinel,
                    lon_sentinel,
                    elev_nan,
                    mat_nan,
                    if (row_prev_iy != null) &row_prev else null,
                    &row_curr,
                    &row_next,
                    row_curr_iy.?,
                    self.lat_ud,
                    self.lon_ud,
                    self.elev,
                    self.matc
                );

                const tmp = row_prev;
                row_prev = row_curr;
                row_curr = row_next;
                row_next = tmp;
                row_prev_iy = row_curr_iy;
                row_curr_iy = row_next_iy;
                row_next_iy = null;
                row_next.reset();
            }
        }

        if (row_curr_iy) |iy| {
            try Helpers.processRow(
                handler,
                tile_stride,
                halo_stride,
                tiles_x,
                lat_sentinel,
                lon_sentinel,
                elev_nan,
                mat_nan,
                if (row_prev_iy != null) &row_prev else null,
                &row_curr,
                null,
                iy,
                self.lat_ud,
                self.lon_ud,
                self.elev,
                self.matc
            );
        }
    }
};
