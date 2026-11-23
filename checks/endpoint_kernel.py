#!/usr/bin/env python3
"""Calibrate the root and frequency identities used for the stopped collar."""

import math


def report_variance(noise):
    points = range(3)
    positions = (-1.0, 0.0, 1.0)

    def root(x, y):
        return 1.0 + 0.6 * positions[x] * positions[y]

    def channel(visible, hidden):
        return (1.0 - noise) if visible == hidden else (noise / 2.0)

    variance, root_error = 0.0, 0.0
    for u in points:
        for v in points:
            weights = [(channel(u, x) * channel(v, y), root(x, y))
                       for x in points for y in points]
            mean = sum(w * a for w, a in weights)
            square = sum(w * a * a for w, a in weights)
            variance += (square - mean * mean) / 9.0
            root_error += sum(w * (a - math.sqrt(square)) ** 2
                              for w, a in weights) / 9.0
    assert root_error <= 2.0 * variance + 1e-14
    return variance, root_error


def symbol_model(q, twist):
    variance = 0.0
    logdet = 0.0
    for j in range(q):
        angle = (2.0 * math.pi * j + twist) / q
        laplace = 4.0 * math.sin(angle / 2.0) ** 2
        symbol = 1.0 + 0.4 * laplace
        root = 1.0 / math.sqrt(symbol)
        variance += (root - 1.0) ** 2 / (laplace + q ** -2)
        logdet += math.log(root)
    return variance / (q * q), logdet


for noise in [0.2, 0.05, 0.01, 0.001]:
    variance, root_error = report_variance(noise)
    print(f"noise={noise:.3g}  variance={variance:.10g}  "
          f"squared root error={root_error:.10g}  ratio={root_error / variance:.7g}")

for q in [32, 128, 512, 2048]:
    variance, determinant = symbol_model(q, 0.0)
    twisted, twisted_determinant = symbol_model(q, 2.1)
    print(f"q={q:4d}  q×variance={q * variance:.12g}  "
          f"twisted={q * twisted:.12g}  determinant change="
          f"{abs(determinant - twisted_determinant):.3g}")
