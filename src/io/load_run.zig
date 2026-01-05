const std = @import("std");
const geo_attr = @import("geo_attr");
const utils = @import("utils");

pub const LoadRun = struct {
    runfile: []const u8 = undefined,
    buf_reader: *std.Io.Reader = undefined,
    err_log: *std.Io.Writer = undefined,

    pub fn load(self: *LoadRun) !void {
        var lat_lon_rng_n_tile_specs = geo_attr.LatLonRangeAndTileSpecs{ .reader = self.buf_reader, .err_log = self.err_log };
        try lat_lon_rng_n_tile_specs.load(self.runfile);
        var sim_init = geo_attr.SimInit{ .reader = self.buf_reader, .err_log = self.err_log };
        try sim_init.load(self.runfile);
        var bin_writer = utils.FileWriter{ .err_log = self.err_log, .is_err_log = false };
        try bin_writer.create("tiles/geo_attr.bin");
        defer bin_writer.close();
        bin_writer.writer();
        var geo_attributes = geo_attr.GeoAttr{ .reader = self.buf_reader, .err_log = self.err_log, .bin_writer = bin_writer.buf_writer };
        try geo_attributes.loadAndWriteBin(&lat_lon_rng_n_tile_specs, self.runfile);
        try geo_attributes.bin_writer.flush();
        var geo_attr_bin_data = geo_attr.GeoAttrBinData{};
        try geo_attr_bin_data.readBinToArraysWithArena(&geo_attributes, "tiles/geo_attr.bin");
        defer geo_attr_bin_data.deinit();
        std.debug.print("test tilespec print: lat_max_ud: {d}, simulation start year: {d}, lat: {d}, lon: {d}, elev: {d}, mat: {d}, ix: {d}, iy: {d}, nx: {d}, ny: {d}\n", .{ lat_lon_rng_n_tile_specs.lat_max_ud, sim_init.start_yr, geo_attributes.lat_ud, geo_attributes.lon_ud, geo_attributes.elevation, geo_attributes.matc, geo_attributes.ix, geo_attributes.iy, geo_attributes.nx, geo_attributes.ny });
    }
};
