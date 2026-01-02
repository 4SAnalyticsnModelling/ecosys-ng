const std = @import("std");
const utils = @import("utils");
const parser = @import("input_parser");
const max_path_len = 1024;
///Lat-lon ranges and tile sepcifications
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
    }
};
///Geographical attributes
pub const GeoAttr = packed struct {
    reader: *std.Io.Reader = undefined,
    err_log: *std.Io.Writer = undefined,
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

    pub fn load(self: *GeoAttr, lat_lon_rng_n_tile_specs: *const LatLonRangeAndTileSpecs, file_name: []const u8) !void {
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
};
