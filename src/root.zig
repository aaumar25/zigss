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
            const magnitude = std.math.sqrt(
                std.math.pow(f32, a.magnitude, 2) +
                    std.math.pow(f32, b.magnitude, 2),
            );
            const phase = std.math.atan2(-b.magnitude, a.magnitude);
            return .{
                .frequency = a.frequency,
                .magnitude = magnitude,
                .phase = phase,
            };
        } else return error.NotSupported;
    }
};

pub const Complex = union(u1) {
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

        pub fn toPolar(self: Cartesian) Polar {
            return .{
                .r = @sqrt(self.re * self.re + self.im * self.im),
                .theta = std.math.atan2(self.re, self.im),
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

        pub fn toCartesian(self: Polar) Cartesian {
            return .{
                .a = self.r * std.math.cos(self.theta),
                .b = self.r * std.math.sin(self.theta),
            };
        }
    };

    /// Calculate the reciprocal of a complex number
    pub fn reciprocal(input: Complex) Complex {
        switch (input) {
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
