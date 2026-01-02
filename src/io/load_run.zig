const std = @import("std");
const geo_attr = @import("geo_attr");

pub const LoadRun = struct {
    runfile: []const u8 = undefined,
    buf_reader: *std.Io.Reader = undefined,
    err_log: *std.Io.Writer = undefined,

    pub fn load(self: *LoadRun) !void {
        var lat_lon_rng_n_tile_specs = geo_attr.LatLonRangeAndTileSpecs{ .reader = self.buf_reader, .err_log = self.err_log };
        try lat_lon_rng_n_tile_specs.load(self.runfile);
        var sim_init = geo_attr.SimInit{ .reader = self.buf_reader, .err_log = self.err_log };
        try sim_init.load(self.runfile);
        var geo_attributes = geo_attr.GeoAttr{ .reader = self.buf_reader, .err_log = self.err_log };
        try geo_attributes.load(&lat_lon_rng_n_tile_specs, self.runfile);
        std.debug.print("test tilespec print: lat_max_ud: {d}, simulation start year: {d}, lat: {d}, lon: {d}, elev: {d}, mat: {d}, ix: {d}, iy: {d}, nx: {d}, ny: {d}\n", .{ lat_lon_rng_n_tile_specs.lat_max_ud, sim_init.start_yr, geo_attributes.lat_ud, geo_attributes.lon_ud, geo_attributes.elevation, geo_attributes.matc, geo_attributes.ix, geo_attributes.iy, geo_attributes.nx, geo_attributes.ny });
    }
};
