const std = @import("std");

/// SOIL.F writes WOUTS/WOUTP/WOUTQ after the final hourly cycle of a day when
/// `(I/KOUT)*KOUT == I`, and unconditionally at the final day of each scene.
/// The checkpoint option still gates both cases and KOUT=0 disables output.
pub fn shouldPublish(enabled: bool, interval_days: u16, day_of_year: u16, completed_hour: u8, is_scene_final_day: bool) !bool {
    if (day_of_year == 0 or day_of_year > 366 or completed_hour > 24) return error.InvalidCheckpointScheduleInstant;
    if (!enabled or interval_days == 0) return false;
    // Native modern streams may use 0..23; supplied Fortran weather uses
    // source clock hours 1..24. Both represent the same completed day here.
    if (completed_hour != 23 and completed_hour != 24) return false;
    return is_scene_final_day or day_of_year % interval_days == 0;
}

/// A scene requests ROUTS/ROUTP only at its initial boundary. The executable
/// must restore the complete committed bundle before any management or hourly
/// science is applied.
pub fn shouldRestore(resume_enabled: bool, scene_hour_index: usize) bool {
    return resume_enabled and scene_hour_index == 0;
}

test "checkpoint cadence preserves SOIL daily modulo and final-scene rules" {
    try std.testing.expect(!try shouldPublish(true, 10, 20, 22, false));
    try std.testing.expect(try shouldPublish(true, 10, 20, 23, false));
    try std.testing.expect(try shouldPublish(true, 10, 20, 24, false));
    try std.testing.expect(!try shouldPublish(true, 10, 21, 23, false));
    try std.testing.expect(try shouldPublish(true, 10, 21, 23, true));
    try std.testing.expect(!try shouldPublish(false, 10, 20, 23, true));
    try std.testing.expect(!try shouldPublish(true, 0, 20, 23, true));
}

test "checkpoint restore occurs only before the first scene hour" {
    try std.testing.expect(shouldRestore(true, 0));
    try std.testing.expect(!shouldRestore(true, 1));
    try std.testing.expect(!shouldRestore(false, 0));
}
