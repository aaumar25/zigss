const std = @import("std");
const zigss = @import("zigss.zig");

/// Since the length of the wave may be known in runtime, all functions in wave
/// return an n-dimensional slice of f32.
pub const wave = struct {
    /// Given an n-dimensional of  slice, array, or vector of f32, return an
    /// n-dimensional slice of f32. If given an f32, return a 1-dimensional
    /// slice of f32. The caller own the memory, and required to release the
    /// memory.
    pub fn sine(
        input: anytype,
        allocator: std.mem.Allocator,
    ) ![]f32 {
        const ti = @typeInfo(@TypeOf(input));
        switch (ti) {
            inline .float,
            .comptime_float,
            => {},
            inline .pointer => |ptr| {
                const result = try allocator.dupe(
                    ptr.child,
                    input,
                );
                for (result) |*res| {
                    res.* = @sin(res.*);
                }
                return result;
            },
            inline .array => {
                const result = try allocator.dupe(
                    f32,
                    &input,
                );
                for (result) |*res| {
                    res.* = @sin(res.*);
                }
                return result;
            },
            inline .vector => {},
            inline else => return error.InvalidInput,
        }
    }

    test sine {
        // From array
        const allocator = std.testing.allocator;
        const radians = [5]f32{
            0.0,
            zigss.pi / 6.0,
            zigss.pi / 4.0,
            zigss.pi / 3.0,
            zigss.pi / 2.0,
        };
        const expect = [5]f32{ 0.0, 0.5, 0.707, 0.866, 1 };
        {
            const result = try wave.sine(radians, allocator);
            defer allocator.free(result);

            for (result, expect) |r, e| {
                try std.testing.expectApproxEqAbs(e, r, 0.001);
            }
        }
        // From slice
        {
            var slice = try allocator.dupe(f32, &radians);
            defer allocator.free(slice);
            const result = try wave.sine(&slice, allocator);
            defer allocator.free(result);
            for (result, expect) |r, e| {
                try std.testing.expectApproxEqAbs(e, r, 0.001);
            }
        }
    }
};
