//! This library is built using the following books:
//! [1] B. P. Lathi and R. A. Green, Linear Systems and Signals. New York:
//!     Oxford University Press, 2018.
const std = @import("std");
const pi: f32 = @floatCast(std.math.pi);

/// Sinusoids can be formulated with radian frequency or hertz frequency.
/// For simplicity, this library will represent the sinusoids in radian.
/// The sinusoids is represented as
///     x(t) = C cos (wt + θ),
/// where:
///     - C: magnitude of sinusoids
///     - w: radian frequency, equal to 2πf, with f as hertz frequency
///     - θ: phase (radian)
///
/// Construct the sinusoids signal using initRadianFreq() or initHzFreq()
pub const Sinusoid = struct {
    magnitude: f32,
    frequency: f32,
    phase: f32,

    pub fn initRadianFreq(magnitude: f32, frequency: f32, phase: f32) Sinusoid {
        return .{
            .magnitude = magnitude,
            .frequency = frequency,
            .phase = phase,
        };
    }

    pub fn initHzFreq(magnitude: f32, frequency: f32, phase: f32) Sinusoid {
        return .{
            .magnitude = magnitude,
            .frequency = 2.0 * pi * frequency,
            .phase = phase,
        };
    }

    pub fn add(a: Sinusoid, b: Sinusoid) !Sinusoid {
        if (a.frequency == b.frequency) {
            const magnitude = @sqrt(a.magnitude * a.magnitude +
                b.magnitude * b.magnitude);
            const phase = std.math.atan2(-b.magnitude, a.magnitude);
            return .{
                .frequency = a.frequency,
                .magnitude = magnitude,
                .phase = phase,
            };
        } else return error.NotSupported;
    }
};

pub const Complex = union(enum) {
    cartesian: Cartesian,
    polar: Polar,
    /// Complex number z can be represented in Cartesian form as follows:
    ///         z = a + jb,
    /// where a and b (the abscissa and the ordinate) of z are the real part
    /// and the imaginary part, respectively, of z.
    pub const Cartesian = struct {
        /// Real part of complex numbers
        re: f32,
        /// Imaginary part of complex numbers
        im: f32,

        pub fn toPolar(z: Cartesian) Polar {
            return .{
                .r = @sqrt(z.re * z.re + z.im * z.im),
                .theta = std.math.atan2(z.im, z.re),
            };
        }
    };

    /// Complex numbers may also be expressed in terms of polar coordinates.
    /// if (r, θ) are the polar coordinates of a point z = a + jb, then
    ///     a = r cos(θ)    and     b = r sin(θ)
    /// Consequently,
    ///     z = r(cos(θ) + j sin(θ)) = r e^(j θ) by Euler's formula
    pub const Polar = struct {
        /// Distance of complex numbers to origin
        r: f32,
        /// Angle of complex numbers
        theta: f32,

        pub fn toCartesian(z: Polar) Cartesian {
            return .{
                .re = z.r * @cos(z.theta),
                .im = z.r * @sin(z.theta),
            };
        }
    };

    /// Calculate the reciprocal of a complex number
    pub fn reciprocal(z: Complex) Complex {
        switch (z) {
            .cartesian => |cartesian| {
                const polar = cartesian.toPolar();
                return .{
                    .polar = .{ .r = 1.0 / polar.r, .theta = -polar.theta },
                };
            },
            .polar => |polar| return .{
                .polar = .{ .r = 1.0 / polar.r, .theta = -polar.theta },
            },
        }
    }

    /// Calculate the conjugate of a complex number
    pub fn conjugate(z: Complex) Complex {
        switch (z) {
            .cartesian => |cartesian| return .{
                .cartesian = .{ .re = cartesian.re, .im = -cartesian.im },
            },
            .polar => |polar| return .{
                .polar = .{ .r = polar.r, .theta = -polar.theta },
            },
        }
    }
};

pub fn degreeToRadian(degree: f32) f32 {
    return degree * pi / 180.0;
}

pub fn radianToDegree(radian: f32) f32 {
    return radian * 180.0 / pi;
}

test "addSinusoid" {
    const a: Sinusoid = .initRadianFreq(1, pi, 0);
    const b: Sinusoid = .initRadianFreq(-std.math.sqrt(3.0), pi, 0);
    const expected = try Sinusoid.add(a, b);
    const actual: Sinusoid = .initRadianFreq(2, pi, degreeToRadian(60));
    try std.testing.expectEqual(expected.magnitude, actual.magnitude);
    try std.testing.expectEqual(expected.frequency, actual.frequency);
    try std.testing.expectEqual(expected.phase, actual.phase);
}

test "ConvertComplexNumbers" {
    const a: Complex = .{ .cartesian = .{ .re = 2, .im = 3 } };
    const b: Complex = .{ .cartesian = .{ .re = -2, .im = 1 } };
    const c: Complex = .{ .cartesian = .{ .re = -2, .im = -3 } };
    const d: Complex = .{ .cartesian = .{ .re = 1, .im = -3 } };
    const expected_a: Complex = .{
        .polar = .{
            .r = @sqrt(13.0),
            .theta = std.math.atan2(@as(f32, 3.0), @as(f32, 2.0)),
        },
    };
    const expected_b: Complex = .{
        .polar = .{
            .r = @sqrt(5.0),
            .theta = std.math.atan2(@as(f32, 1.0), -@as(f32, 2.0)),
        },
    };
    const expected_c: Complex = .{
        .polar = .{
            .r = @sqrt(13.0),
            .theta = std.math.atan2(-@as(f32, 3.0), -@as(f32, 2.0)),
        },
    };
    const expected_d: Complex = .{
        .polar = .{
            .r = @sqrt(10.0),
            .theta = std.math.atan2(-@as(f32, 3.0), @as(f32, 1.0)),
        },
    };
    // Testing cartesian to polar
    try std.testing.expectEqual(a.cartesian.toPolar(), expected_a.polar);
    try std.testing.expectEqual(b.cartesian.toPolar(), expected_b.polar);
    try std.testing.expectEqual(c.cartesian.toPolar(), expected_c.polar);
    try std.testing.expectEqual(d.cartesian.toPolar(), expected_d.polar);
    // Testing polar to cartesian

    try structExpectAproxEqAbs(a.cartesian, expected_a.polar.toCartesian(), 1e-6);
    try structExpectAproxEqAbs(b.cartesian, expected_b.polar.toCartesian(), 1e-6);
    try structExpectAproxEqAbs(c.cartesian, expected_c.polar.toCartesian(), 1e-6);
    try structExpectAproxEqAbs(d.cartesian, expected_d.polar.toCartesian(), 1e-6);
}

fn structExpectAproxEqAbs(expected: anytype, actual: anytype, tolerance: anytype) !void {
    const ti = @typeInfo(@TypeOf(expected));
    if (ti != .@"struct") @compileError("UnsupportedType");
    inline for (ti.@"struct".fields) |field| {
        const e_val = @field(expected, field.name);
        const a_val = @field(actual, field.name);
        try std.testing.expectApproxEqAbs(e_val, a_val, tolerance);
    }
}
