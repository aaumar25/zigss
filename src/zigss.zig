//! Root of zigss library
const std = @import("std");
pub const signal = @import("signal.zig");

pub const pi = 3.1415927;

const _ = signal.wave.sine(0.1, std.testing.allocator);

test {
    std.testing.refAllDeclsRecursive(@This());
}
