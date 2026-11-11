"""Shortest generators of the central flux classes and their plaquette costs.

The class-k element is (e^{2 pi i Y})^k on the five-level code. Its eigenphases, in units of pi, are
-2k/3 on the three colour levels and k on the two doublet levels, each modulo 2. A generator H with
exp(H) equal to that element and trace zero has eigenvalues i pi x_j with x_j a lift of the phase. The
shortest one minimises tr_5 H^dag H = pi^2 sum x_j^2, and the coefficient of the transmutation exponent
for alpha = g^2 / 4 pi is beta_k = 8 pi^2 / tr_5 H^dag H = 8 / sum x_j^2.
"""
from fractions import Fraction as F
from itertools import product

SPAN = 3


def class_phases(k):
    return [F(-2 * k, 3)] * 3 + [F(k)] * 2


def shortest_generator(k):
    base = class_phases(k)
    best = None
    for shifts in product(range(-SPAN, SPAN + 1), repeat=5):
        x = [b + 2 * n for b, n in zip(base, shifts)]
        if sum(x) != 0:
            continue
        s = sum(v * v for v in x)
        if best is None or s < best[0] or (s == best[0] and x < best[1]):
            best = (s, x)
    return best


def costs():
    rows = {}
    for k in range(6):
        s, x = shortest_generator(k)
        rows[str(k)] = dict(eigenphases_over_pi=[str(v) for v in x], trace_over_pi_squared=str(s),
                            beta=None if s == 0 else str(F(8) / s))
    return rows


def checks(rows):
    beta = {k: (None if r['beta'] is None else F(r['beta'])) for k, r in rows.items()}
    weak = [F(v) for v in rows['3']['eigenphases_over_pi']]
    s3 = F(rows['3']['trace_over_pi_squared'])
    # exponent identity 2 pi / (beta alpha) = tr H^dag H / g^2 at alpha = g^2 / (4 pi), in units pi^2 / g^2
    exponent_identity = all(F(8) / beta[k] == F(rows[k]['trace_over_pi_squared']) for k in '12345')
    return dict(
        betas_are_12_5_3_4_3_12_5=[beta[k] for k in '12345'] == [F(12, 5), 3, 4, 3, F(12, 5)],
        class_zero_has_zero_generator=rows['0']['trace_over_pi_squared'] == '0',
        weak_centre_generator_is_plus_or_minus_2_pi_i_T3=weak in ([0, 0, 0, 1, -1], [0, 0, 0, -1, 1]),
        weak_centre_trace_is_2_pi_squared=s3 == 2,
        weak_centre_element_is_the_2_pi_rotation=[(-1) ** int(v) for v in weak] == [1, 1, 1, -1, -1],
        weak_centre_beta_equals_N_c_plus_one=beta['3'] == 3 + 1,
        exponent_equals_gaussian_plaquette_action=exponent_identity,
        conjugate_classes_have_equal_costs=all(beta[str(k)] == beta[str(6 - k)] for k in (1, 2)),
        lifts_inside_search_window=all(abs(F(v)) < 2 * SPAN - 1 for r in rows.values()
                                       for v in r['eigenphases_over_pi']))
