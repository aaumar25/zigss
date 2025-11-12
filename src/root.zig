//! This library is built using the following books:
//! [1] B. P. Lathi and R. A. Green, Linear Systems and Signals. New York:
//!     Oxford University Press, 2018.
const std = @import("std");
const pi: f64 = @floatCast(std.math.pi);

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
    magnitude: f64,
    frequency: f64,
    phase: f64,

    pub fn initRadianFreq(magnitude: f64, frequency: f64, phase: f64) Sinusoid {
        return .{
            .magnitude = magnitude,
            .frequency = frequency,
            .phase = phase,
        };
    }

    pub fn initHzFreq(magnitude: f64, frequency: f64, phase: f64) Sinusoid {
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
    // TODO: Can we represent the complex number in a function? See Example
    // B.5.
    cartesian: Cartesian,
    polar: Polar,
    /// Complex number z can be represented in Cartesian form as follows:
    ///         z = a + jb,
    /// where a and b (the abscissa and the ordinate) of z are the real part
    /// and the imaginary part, respectively, of z.
    pub const Cartesian = struct {
        /// Real part of complex numbers
        re: f64,
        /// Imaginary part of complex numbers
        im: f64,

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
        r: f64,
        /// Angle of complex numbers
        theta: f64,

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

    /// Addition of complex numbers. Return the complex number in cartesian
    /// form.
    pub fn add(a: Complex, b: Complex) Complex {
        // Ensure that a and b are in cartesian for ease of calculation.
        const a_cart: Complex = .{
            .cartesian = switch (a) {
                .cartesian => |cartesian| cartesian,
                .polar => |polar| polar.toCartesian(),
            },
        };
        const b_cart: Complex = .{
            .cartesian = switch (b) {
                .cartesian => |cartesian| cartesian,
                .polar => |polar| polar.toCartesian(),
            },
        };
        return .{
            .cartesian = .{
                .re = a_cart.cartesian.re + b_cart.cartesian.re,
                .im = a_cart.cartesian.im + b_cart.cartesian.im,
            },
        };
    }

    /// Substraction of complex numbers. Return the complex number in
    /// cartesian form.
    pub fn sub(a: Complex, b: Complex) Complex {
        // Ensure that a and b are in cartesian for ease of calculation.
        const a_cart: Complex = .{
            .cartesian = switch (a) {
                .cartesian => |cartesian| cartesian,
                .polar => |polar| polar.toCartesian(),
            },
        };
        const b_cart: Complex = .{
            .cartesian = switch (b) {
                .cartesian => |cartesian| cartesian,
                .polar => |polar| polar.toCartesian(),
            },
        };
        return .{
            .cartesian = .{
                .re = a_cart.cartesian.re - b_cart.cartesian.re,
                .im = a_cart.cartesian.im - b_cart.cartesian.im,
            },
        };
    }

    /// Multiplication of complex numbers (a * b). Return the complex number
    /// in cartesian form. Support multiplication with scalar or other complex
    /// number.
    pub fn mul(
        a: Complex,
        b: union(enum) { complex: Complex, scalar: f64 },
    ) Complex {
        switch (b) {
            .complex => |b_val| {
                // Ensure that a and b are in polar for ease of calculation.
                const a_p: Complex = .{
                    .polar = switch (a) {
                        .cartesian => |cartesian| cartesian.toPolar(),
                        .polar => |polar| polar,
                    },
                };
                const b_p: Complex = .{
                    .polar = switch (b_val) {
                        .cartesian => |cartesian| cartesian.toPolar(),
                        .polar => |polar| polar,
                    },
                };
                const result_polar: Complex.Polar = .{
                    .r = a_p.polar.r * b_p.polar.r,
                    .theta = a_p.polar.theta + b_p.polar.theta,
                };
                return .{ .cartesian = result_polar.toCartesian() };
            },
            .scalar => |b_val| {
                switch (a) {
                    .cartesian => |a_val| {
                        return .{ .cartesian = .{
                            .re = a_val.re * b_val,
                            .im = a_val.im * b_val,
                        } };
                    },
                    .polar => |a_val| {
                        const a_cart = a_val.toCartesian();
                        return .{ .cartesian = .{
                            .re = a_cart.re * b_val,
                            .im = a_cart.im * b_val,
                        } };
                    },
                }
            },
        }
    }

    /// Division of complex numbers (a / b). Return the complex number in
    /// cartesian form.
    pub fn div(a: Complex, b: Complex) Complex {
        // Ensure that a and b are in polar for ease of calculation.
        const a_p: Complex = .{
            .polar = switch (a) {
                .cartesian => |cartesian| cartesian.toPolar(),
                .polar => |polar| polar,
            },
        };
        const b_p: Complex = .{
            .polar = switch (b) {
                .cartesian => |cartesian| cartesian.toPolar(),
                .polar => |polar| polar,
            },
        };
        const result_polar: Complex.Polar = .{
            .r = a_p.polar.r / b_p.polar.r,
            .theta = a_p.polar.theta - b_p.polar.theta,
        };
        return .{ .cartesian = result_polar.toCartesian() };
    }

    /// Powers of complex number a with b (a ^ b).
    pub fn pow(a: Complex, b: f64) Complex {
        // Ensure that a is in polar for ease of calculation.
        const a_p: Complex = .{
            .polar = switch (a) {
                .cartesian => |cartesian| cartesian.toPolar(),
                .polar => |polar| polar,
            },
        };
        const result_polar: Complex.Polar = .{
            .r = std.math.pow(f64, a_p.polar.r, b),
            .theta = a_p.polar.theta * b,
        };
        return .{ .cartesian = result_polar.toCartesian() };
    }

    /// Roots of complex number a with b (a ^ (1 / b)).
    pub fn roots(a: Complex, b: f64) Complex {
        // NOTE: This function only return principal value of root complex.
        // If it is required to return all the rots, see B.12 formula.
        // Ensure that a is in polar for ease of calculation.
        const a_p: Complex = .{
            .polar = switch (a) {
                .cartesian => |cartesian| cartesian.toPolar(),
                .polar => |polar| polar,
            },
        };
        const result_polar: Complex.Polar = .{
            .r = std.math.pow(f64, a_p.polar.r, 1 / b),
            .theta = a_p.polar.theta / b,
        };
        return .{ .cartesian = result_polar.toCartesian() };
    }

    /// Natural logarithm of a complex number. Return the complex number in
    /// cartesian form.
    pub fn log(a: Complex) Complex {
        // Ensure that a is in polar for ease of calculation.
        const a_p: Complex = .{
            .polar = switch (a) {
                .cartesian => |cartesian| cartesian.toPolar(),
                .polar => |polar| polar,
            },
        };
        const result_polar: Complex.Polar = .{
            .r = @log(a_p.polar.r),
            .theta = a_p.polar.theta,
        };
        return .{ .cartesian = result_polar.toCartesian() };
    }
};

pub fn degreeToRadian(degree: f64) f64 {
    return degree * pi / 180.0;
}

pub fn radianToDegree(radian: f64) f64 {
    return radian * 180.0 / pi;
}

test "addSinusoid" {
    const a: Sinusoid = .initRadianFreq(1, pi, 0);
    const b: Sinusoid = .initRadianFreq(-std.math.sqrt(3.0), pi, 0);
    const expected = try Sinusoid.add(a, b);
    const actual: Sinusoid = .initRadianFreq(2, pi, degreeToRadian(60));
    try std.testing.expectApproxEqAbs(expected.magnitude, actual.magnitude, 1e-6);
    try std.testing.expectApproxEqAbs(expected.frequency, actual.frequency, 1e-6);
    try std.testing.expectApproxEqAbs(expected.phase, actual.phase, 1e-6);
}

test "Convert Complex Numbers" {
    const a: Complex = .{ .cartesian = .{ .re = 2, .im = 3 } };
    const b: Complex = .{ .cartesian = .{ .re = -2, .im = 1 } };
    const c: Complex = .{ .cartesian = .{ .re = -2, .im = -3 } };
    const d: Complex = .{ .cartesian = .{ .re = 1, .im = -3 } };
    const expected_a: Complex = .{
        .polar = .{
            .r = @sqrt(13.0),
            .theta = std.math.atan2(@as(f64, 3.0), @as(f64, 2.0)),
        },
    };
    const expected_b: Complex = .{
        .polar = .{
            .r = @sqrt(5.0),
            .theta = std.math.atan2(@as(f64, 1.0), -@as(f64, 2.0)),
        },
    };
    const expected_c: Complex = .{
        .polar = .{
            .r = @sqrt(13.0),
            .theta = std.math.atan2(-@as(f64, 3.0), -@as(f64, 2.0)),
        },
    };
    const expected_d: Complex = .{
        .polar = .{
            .r = @sqrt(10.0),
            .theta = std.math.atan2(-@as(f64, 3.0), @as(f64, 1.0)),
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

test "Arithmatical Operation of Complex Numbers" {
    const z1: Complex = .{ .cartesian = .{ .re = 3, .im = 4 } };
    const z2: Complex = .{ .cartesian = .{ .re = 2, .im = 3 } };
    // Multiplication
    {
        const actual = z1.mul(.{ .complex = z2 }).cartesian;
        try structExpectAproxEqAbs(
            Complex.Cartesian{ .re = -6, .im = 17 },
            actual,
            1e-6,
        );
    }
    // Division
    try structExpectAproxEqAbs(
        Complex.Cartesian{ .re = 18.0 / 13.0, .im = -1.0 / 13.0 },
        z1.div(z2).cartesian,
        1e-6,
    );
    // Addition
    try structExpectAproxEqAbs(
        Complex.Cartesian{ .re = 5, .im = 7 },
        z1.add(z2).cartesian,
        1e-6,
    );
    // Substraction
    try structExpectAproxEqAbs(
        Complex.Cartesian{ .re = 1, .im = 1 },
        z1.sub(z2).cartesian,
        1e-6,
    );
    // 2nd test
    const z3: Complex = .{ .polar = .{ .r = 2, .theta = pi / 4.0 } };
    const z4: Complex = .{ .polar = .{ .r = 8, .theta = pi / 3.0 } };
    // Substraction
    {
        const actual = z3.mul(.{ .scalar = 2.0 }).sub(z4);
        try structExpectAproxEqAbs(
            Complex.Cartesian{ .re = 2.0 * @sqrt(2.0) - 4.0, .im = 2.0 * @sqrt(2.0) - 4.0 * @sqrt(3.0) },
            actual.cartesian,
            1e-6,
        );
    }
    // Division and Power
    {
        const actual = z3.div(z4.pow(2.0));
        const expected: Complex = .{
            .polar = .{ .r = 1.0 / 32.0, .theta = -5.0 * pi / 12.0 },
        };
        try structExpectAproxEqAbs(
            expected.polar.toCartesian(),
            actual.cartesian,
            1e-6,
        );
    }
    // Roots
    {
        const actual = z4.roots(3.0);
        const expected: Complex = .{ .polar = .{ .r = 2, .theta = pi / 9.0 } };
        try structExpectAproxEqAbs(
            expected.polar.toCartesian(),
            actual.cartesian,
            1e-6,
        );
    }
}

/// This function is intended for testing only. It checks every value in a
/// struct whether the expected value is equal to an actual value with
/// absoulte tolerance.
fn structExpectAproxEqAbs(
    expected: anytype,
    actual: anytype,
    tolerance: anytype,
) !void {
    const ti = @typeInfo(@TypeOf(expected));
    if (ti != .@"struct") @compileError("UnsupportedType");
    inline for (ti.@"struct".fields) |field| {
        const e_val = @field(expected, field.name);
        const a_val = @field(actual, field.name);
        try std.testing.expectApproxEqAbs(e_val, a_val, tolerance);
    }
}
