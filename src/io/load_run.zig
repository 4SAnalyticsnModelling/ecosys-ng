const std = @import("std");
const geo_attr = @import("geo_attr");
const utils = @import("utils");

pub const LoadRun = struct {
    allocator: std.mem.Allocator = undefined,
    runfile: []const u8 = undefined,
    buf_reader: *std.Io.Reader = undefined,
    err_log: *std.Io.Writer = undefined,
    geo_attr_bin_path: []const u8 = "tiles/geo_attr.bin",

    pub fn load(self: *LoadRun) !void {
        var lat_lon_rng_n_tile_specs = geo_attr.LatLonRangeAndTileSpecs{
            .file_name = self.runfile,
            .reader = self.buf_reader,
            .err_log = self.err_log,
        };
        try lat_lon_rng_n_tile_specs.load();
        var sim_init = geo_attr.SimInit{
            .file_name = self.runfile,
            .reader = self.buf_reader,
            .err_log = self.err_log,
        };
        try sim_init.load();
        var bin_writer = utils.FileWriter{
            .err_log = self.err_log,
            .is_err_log = false,
        };
        try bin_writer.create(self.geo_attr_bin_path);
        defer bin_writer.close();
        bin_writer.writer();
        var geo_attributes = geo_attr.GeoAttr{
            .file_name = self.runfile,
            .allocator = self.allocator,
            .reader = self.buf_reader,
            .err_log = self.err_log,
            .bin_writer = bin_writer.buf_writer,
        };
        try geo_attributes.loadAndWriteBin(&lat_lon_rng_n_tile_specs);
        try geo_attributes.bin_writer.flush();
        var geo_attr_bin_data = geo_attr.GeoAttrBinData{
            .allocator = self.allocator,
            .err_log = self.err_log,
        };
        try geo_attr_bin_data.readBinByTileAndFillArrays(&geo_attributes, &lat_lon_rng_n_tile_specs, self.geo_attr_bin_path, LoadRun.handleGeoAttrTile);
        var tile_reader = try geo_attr.GeoAttrBinData.TileReader.init(self.allocator, self.err_log, &lat_lon_rng_n_tile_specs, self.geo_attr_bin_path);
        defer tile_reader.deinit();
        while (true) {
            const batch_opt = try tile_reader.next();
            if (batch_opt == null) break;
            const batch = batch_opt.?;
            std.debug.print("tile ({d}, {d}) count={d}\n", .{ batch.tile_ix, batch.tile_iy, batch.records.len });
        }
    }

    fn handleGeoAttrTile(tile_ix: usize, tile_iy: usize, tile_nx: usize, tile_ny: usize, halo_stride: usize, tile_lat_ud: []const i32, tile_lon_ud: []const i32, tile_elev: []const f32, tile_matc: []const f32, records: []const geo_attr.GeoAttr.Record) !void {
        _ = records;
        std.debug.print(
            "test tile ({d}, {d}) tile_nx={d} tile_ny={d} stride={d}\n",
            .{ tile_ix, tile_iy, tile_nx, tile_ny, halo_stride },
        );
        var y: usize = 0;
        while (y < tile_ny) : (y += 1) {
            var x: usize = 0;
            while (x < tile_nx) : (x += 1) {
                const idx = (x + 1) + (y + 1) * halo_stride;
                std.debug.print(
                    "  cell ({d}, {d}) lat_ud={d} lon_ud={d} elev={} matc={}\n",
                    .{ x, y, tile_lat_ud[idx], tile_lon_ud[idx], tile_elev[idx], tile_matc[idx] },
                );
            }
        }
    }
};
