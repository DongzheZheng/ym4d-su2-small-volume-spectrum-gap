"""Exact rational bookkeeping for the fixed-volume weak-coupling proof."""
from fractions import Fraction

f = Fraction
delta = f(1, 1000)

case = {
    "wide_action_decay_beta": f(11, 10) - (f(1, 3) + 4 * f(1, 24)),
    "forward_root_decay_beta": (f(1, 24) + f(1, 3) - delta) / 2,
    "chart_contraction_decay_g": 1 - 4 * f(1, 12),
    "outer_suppression_growth_g_inverse": f(2, 3) - f(5, 9),
    "static_hard_decay_g": 2 * f(2, 3) - f(2, 3) - 3 * f(1, 12),
    "electric_young_decay_g": -f(3, 4) - f(3, 8) + 2 * f(3, 4),
    "central_haar_decay_g": -f(3, 4) + 2 * f(3, 4),
    "normal_marginal_squared_decay_g": 2 * (f(5, 9) - f(1, 12)) - 2 * delta,
    "soft_arms_root_decay_g": (f(3, 4) - 2 * delta) / 2,
}
assert case["wide_action_decay_beta"] == f(3, 5)
assert case["static_hard_decay_g"] == f(5, 12)
assert case["outer_suppression_growth_g_inverse"] == f(1, 9)
assert case["electric_young_decay_g"] == f(3, 8)
assert case["normal_marginal_squared_decay_g"] == f(17, 18) - 2 * delta
assert case["soft_arms_root_decay_g"] / 2 == case["forward_root_decay_beta"]
assert case["forward_root_decay_beta"] > f(1, 8)

upper_cube = f(39, 4) ** 3 * f(30, 13)
odd_lower_cube = 27 * 86
assert upper_cube == f(68445, 32)
assert upper_cube < odd_lower_cube

for label, value in case.items():
    print(label, value)
print("excited_even_upper_cube", upper_cube)
print("odd_lower_cube", odd_lower_cube)
print("slack_over_beta_one_eighth_before_o1", case["forward_root_decay_beta"] - f(1, 8))
