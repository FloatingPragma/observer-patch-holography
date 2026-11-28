"""Offline check of the strength of the registered cubic energy-lapse family."""

import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
FIXTURE_HASH = "578776c047e327661f496a261986274f5e7e0f8f0498f4729c16d73ae269b2ff"
Q = F(1, 8)
MASS, MASS_U = F("9.1093837139e-31"), F("0.0000000028e-31")
BIG_G, BIG_G_U = F("6.67430e-11"), F("0.00015e-11")
PLANCK, C = F("6.62607015e-34"), 299792458
ENERGY_GEV = F("0.00051099895069")
LENGTH_FM = F(236221477831373991, 312500000000000000)
HC_GEV_FM, TARGET_GEV = F("0.1973269804"), F("6.70883e-39")


def arctan_box(x):
    partial = sum((-1) ** i * x ** (2 * i + 1) / (2 * i + 1) for i in range(36))
    return partial, partial + x**73 / 73


def pi_box():
    a, b = arctan_box(F(1, 5))
    c, d = arctan_box(F(1, 239))
    return 16 * a - 4 * d, 16 * b - 4 * c


def digits(q, roots=0):
    with localcontext() as ctx:
        ctx.prec = 30
        x = Decimal(q.numerator) / Decimal(q.denominator)
        for _ in range(roots):
            x = x.sqrt()
        return format(x, ".14E")


def fixture_energy():
    raw = (HERE / "profile_fixture.json").read_bytes()
    if sha256(raw).hexdigest() != FIXTURE_HASH:
        raise ValueError("unexpected input profile")
    archive = json.loads(raw)
    n, den = archive["side"], archive["denominator"]
    points = list(product(range(n), repeat=3))
    f = dict(zip(points, archive["integer_coordinates"], strict=True))
    assert n == 7 and den == 2**40 and f[(3, 3, 3)] == den
    energy = 0
    for point, value in f.items():
        for axis in range(3):
            end = list(point)
            end[axis] += 1
            energy += (value - f.get(tuple(end), 0)) ** 2
            if point[axis] == 0:
                energy += value**2  # exterior zero on the negative face
    return F(energy, den**2)


def flow_layer(k):
    """Signed-octant average of oriented outward edge currents."""
    count = (k + 1) * (k + 2) // 2
    flow = defaultdict(F)
    octant_energy = F(0)
    for x in range(k + 1):
        for y in range(k - x + 1):
            start = (x, y, k - x - y)
            for axis in range(3):
                end = list(start)
                end[axis] += 1
                weight = F(end[axis], count * (k + 3))
                octant_energy += weight**2
                for signs in product((-1, 1), repeat=3):
                    u = tuple(a * b for a, b in zip(start, signs))
                    v = tuple(a * b for a, b in zip(end, signs))
                    flow[u, v] += weight / 8
    assert octant_energy == F(1, (k + 1) * (k + 3))
    assert sum(value**2 for value in flow.values()) <= octant_energy / 2
    return flow


def check_flow():
    incoming = {(0, 0, 0): F(1)}
    reflected_energy = F(0)
    for k in range(12):
        edges = flow_layer(k)
        outward, next_layer = defaultdict(F), defaultdict(F)
        for (u, v), current in edges.items():
            outward[u] += current
            next_layer[v] += current
            reflected_energy += current**2
        assert outward == incoming and sum(next_layer.values()) == 1
        incoming = next_layer
    # Exact finite-layer sum and analytic tail. Eight reflected octants share
    # each edge with at most three others, including edges on a coordinate axis.
    layers = sum(F(1, (k + 1) * (k + 3)) for k in range(12))
    tail = (F(1, 13) + F(1, 14)) / 2
    assert layers + tail == F(3, 4) and reflected_energy <= layers / 2
    return reflected_energy


def run():
    energy = fixture_energy()
    kappa = 1 / (2 * Q**2 * energy)
    assert kappa == F(4835703278458516698824704, 642754164150366239279661)
    layers = check_flow()
    # Two slopes from a center to grounded endpoints give these exact 1D
    # control energies and exhibit the missing 3D size-independent floor.
    one_dimensional = {
        n: sum((F(1, n)) ** 2 for _ in range(2 * n)) for n in (8, 64, 512)
    }
    assert one_dimensional == {n: F(2, n) for n in one_dimensional}
    pl, ph = pi_box()
    pi = (pl + ph) / 2
    finite_g = Q**2 * energy / (2 * pi)
    g_floor = Q**2 * F(8, 3) / (2 * ph)
    alpha = 2 * pi * BIG_G * MASS**2 / (PLANCK * C)
    alpha_min = 2 * pl * (BIG_G - 5 * BIG_G_U) * (MASS - 5 * MASS_U) ** 2 / (PLANCK * C)
    alpha_max = 2 * ph * (BIG_G + 5 * BIG_G_U) * (MASS + 5 * MASS_U) ** 2 / (PLANCK * C)
    finite_R = finite_g / alpha
    family_R = g_floor / alpha_max
    assert energy >= F(8, 3) and family_R > F("3.78e42")
    assert family_R > F("1.94e21") ** 2 and family_R > F("4.40e10") ** 4
    assert finite_R > F("6.03e42")
    # Suggested changes to q, gap and wave slope enter the invariant in
    # different powers. No common length or energy unit appears here.
    normalized_q_limit_squared = 2 * ph * alpha_max / F(8, 3)
    assert normalized_q_limit_squared < F("6.43e-23") ** 2

    # Check the monomial identity across actual scale changes. The Green
    # coefficient uses the historical GeV/fm records throughout this block.
    def scales(a, u):
        return (
            finite_g * a / (u * HC_GEV_FM * TARGET_GEV),
            u / ENERGY_GEV,
            u * a / HC_GEV_FM,
        )

    reference = finite_g / (TARGET_GEV * ENERGY_GEV**2)
    for a, u in product((LENGTH_FM, 10 * LENGTH_FM), (ENERGY_GEV, 7 * ENERGY_GEV)):
        g, e, v = scales(a, u)
        assert g * e**2 / v == reference
    registered_g, _, registered_light = scales(LENGTH_FM, ENERGY_GEV)
    # Positive control. A system with the matching dimensionless coupling
    # reproduces all three targets when u matches the electron and a matches c.
    matched_g = TARGET_GEV * ENERGY_GEV**2
    matched_length = HC_GEV_FM / ENERGY_GEV
    assert matched_g * matched_length / (ENERGY_GEV * HC_GEV_FM) == TARGET_GEV
    assert ENERGY_GEV * matched_length / HC_GEV_FM == 1
    return {
        "schema": "oph.openwave26.cubic_strength.v1",
        "fixture_sha256": FIXTURE_HASH,
        "dirichlet_energy": str(energy),
        "kappa_exact": str(kappa),
        "q": str(Q),
        "g_fixture": digits(finite_g),
        "reflected_flow_energy_through_layer_11": str(layers),
        "infinite_flow_energy_upper": "3/8",
        "dirichlet_energy_floor": "8/3",
        "g_floor_all_profiles": digits(g_floor),
        "historical_G_pred_over_observed": digits(registered_g),
        "historical_speed_over_c": digits(registered_light),
        "CODATA_electron_gravity_coupling": digits(alpha),
        "CODATA_five_standard_uncertainty_box": [digits(alpha_min), digits(alpha_max)],
        "fixture_dimensionless_mismatch_central": digits(finite_R),
        "family_dimensionless_mismatch_lower": digits(family_R),
        "best_common_error_factor_fixture_electron_fixed": digits(finite_R, 1),
        "best_common_error_factor_fixture_both_units_free": digits(finite_R, 2),
        "best_common_error_factor_family_lower_electron_fixed": digits(family_R, 1),
        "best_common_error_factor_family_lower_both_units_free": digits(family_R, 2),
        "necessary_max_abs_q_times_gap_over_sqrt_speed": digits(
            normalized_q_limit_squared, 1
        ),
        "one_dimensional_control_energies": {
            str(k): str(v) for k, v in one_dimensional.items()
        },
        "common_length_positive_control_fm": digits(matched_length),
    }


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O to enable the scientific checks.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="regenerate the local numerical receipt"
    )
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(path.read_text()):
        raise ArithmeticError("calculated results differ from the archived receipt")
    print(json.dumps(result, indent=2))
