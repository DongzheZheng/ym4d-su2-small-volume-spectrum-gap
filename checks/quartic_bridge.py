#!/usr/bin/env python3
"""Check the scalar OU sample precision and quartic bridge first moment."""

import math
from fractions import Fraction


def ou_precision_ratio(size, mass, twist):
    step = 1.0 / size
    scaled_mass = mass * step
    exact_scale = scaled_mass / math.sinh(scaled_mass)
    exact_mass_square = 4.0 * math.sinh(scaled_mass / 2.0) ** 2
    ratios = []
    for index in range(size):
        angle = (2.0 * math.pi * index + twist) / size
        laplacian = 4.0 * math.sin(angle / 2.0) ** 2
        nearest = laplacian + scaled_mass ** 2
        exact = exact_scale * (laplacian + exact_mass_square)
        ratios.append(exact / nearest)
    return max(abs(value - 1.0) for value in ratios)


def centered_quartic_mean(step, left, right):
    # Boole's exact rule for the two polynomials of degree four.
    points = [Fraction(i, 4) for i in range(5)]
    weights = [Fraction(i, 90) for i in (7, 32, 12, 32, 7)]
    linear = [(1 - s) * left + s * right for s in points]
    quadratic_term = 6 * sum(w * c ** 2 * s * (1 - s)
                             for w, c, s in zip(weights, linear, points))
    quartic_term = 3 * sum(w * s ** 2 * (1 - s) ** 2
                          for w, s in zip(weights, points))
    return step ** 2 * quadratic_term + step ** 3 * quartic_term


for size in (8, 16, 32, 64, 128, 256):
    step = Fraction(1, size)
    twist_ratios = [ou_precision_ratio(size, 0.07, twist)
                    for twist in (0.0, 2.1, math.pi)]
    mean = centered_quartic_mean(step, Fraction(-4, 5), Fraction(11, 10))
    print(f"q={size:3d}  q²×relative precision="
          f"{size ** 2 * max(twist_ratios):.12g}  "
          f"bridge mean={float(mean):.12g}  "
          f"q²×bridge mean={float(mean / step ** 2):.10g}")
