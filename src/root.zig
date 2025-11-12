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
