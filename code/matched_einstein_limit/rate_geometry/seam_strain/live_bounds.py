"""Outward bounds on the primary five-station, volume-matched radar comparison."""

import argparse
from fractions import Fraction
from itertools import combinations
import json

import numpy as np
from mpmath import iv

from ...source_implication.full_algebra.enclosure import Jet
from . import live_radar, live_source, live_survey

KEYS = ("ordinary", "decorrelated", "hidden", "left", "right")


def rational_interval(q):
    return iv.mpf(q.numerator) / q.denominator


def endpoints(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def paired(edges, denominator, state, trace, density_error, effect_error):
    """Enclose relative to the known trace and the fixed timeout range of 3."""
    offsets = []
    for re, im in edges:
        # Subtract 3I before pairing. The true source trace is exact.
        re = re - 3 * denominator * np.eye(8, dtype=int).astype(object)
        radius = Fraction(int(max(np.sum(abs(re) + abs(im), axis=1))), denominator)
        center = (
            sum(
                int(re[i, j]) * Fraction(float(state[j, i].real))
                - int(im[i, j]) * Fraction(float(state[j, i].imag))
                for i in range(8)
                for j in range(8)
            )
            / denominator
            + 3 * trace
        )
        # Every enrolled direction has trace norm <= 2. This includes the
        # commutator with a generator of norm 1/2 and a normalized density.
        error = radius * Fraction(float(density_error)) + 2 * effect_error
        offsets.append(
            rational_interval(center) + iv.mpf([-1, 1]) * rational_interval(error)
        )
    return offsets


def geometry_jets(lengths, variations):
    n = len(KEYS)
    matrix = np.empty((5, 5), dtype=object)
    zero = iv.mpf(0)
    for i in range(5):
        matrix[i, i] = Jet(
            zero, np.full(n, zero, object), np.full((n, n), zero, object)
        )
    for edge, (i, j) in enumerate(combinations(range(5), 2)):
        matrix[i, j] = matrix[j, i] = Jet(
            lengths[edge],
            np.array([variations[k][edge] for k in KEYS], object),
            np.full((n, n), zero, object),
        )
    half = iv.mpf(1) / 2

    def dot(i, j):
        a, b, c = matrix[4, i], matrix[4, j], matrix[i, j]
        return (a * a + b * b - c * c) * half

    area = volume = angle = 0
    for i, j, k in live_survey.FACES:
        a, b = matrix[i, j], matrix[i, k]
        c = matrix[j, k]
        q = (a * a + b * b - c * c) * half
        area += (a * a * b * b - q * q).sqrt() * half
        x, y, z = matrix[4, i], matrix[4, j], matrix[4, k]
        p, q, r = dot(i, j), dot(i, k), dot(j, k)
        det = (
            x * x * y * y * z * z
            + 2 * p * q * r
            - x * x * r * r
            - y * y * q * q
            - z * z * p * p
        )
        root = det.sqrt()
        volume += root * (iv.mpf(1) / 6)
        denom = x * y * z + x * r + y * q + z * p
        if denom.v.a <= 0:
            raise ValueError("solid-angle denominator can vanish")
        reciprocal = denom.unary(1 / denom.v, -1 / (denom.v**2), 2 / (denom.v**3))
        tan = root * reciprocal
        t = tan.v
        angle += 2 * tan.unary(
            iv.atan2(t, 1), 1 / (1 + t * t), -2 * t / (1 + t * t) ** 2
        )
    return area, volume, angle


def charge_rank_at_matched_geometry(fixed):
    """Determinant of the two modular readings on the ordinary and pairing curves."""
    ordinary, pairing = [
        fixed[key]["fixed_outer_area_and_volume"]["modular"]
        for key in ("ordinary", "decorrelated")
    ]
    result = iv.mpf(ordinary["incoming"]) * iv.mpf(pairing["retained"]) - iv.mpf(
        pairing["incoming"]
    ) * iv.mpf(ordinary["retained"])
    assert result.b < 0
    return endpoints(result)


def certificate():
    from . import live_modular

    ports = live_survey.stations()
    q1, den, e1 = live_radar.certified_ranges(1, -1, True, ports)
    q2, den2, e2 = live_radar.certified_ranges(1, 1, True, ports)
    assert den == den2
    kernels = q1 + q2
    effect_error = (e1 + e2) / 2
    source_bounds = live_source.certify()
    material = live_modular.certify()
    errors = source_bounds["channels"]["coupled"]["nominal_state_trace_norm_errors"]
    ref, changes, _ = live_survey.retained(True)
    derr = {k: errors[k] + errors["reference"] for k in ("ordinary", "decorrelated")}
    derr["hidden"] = errors["hidden"]
    derr.update(material["rotation_trace_norm_errors"])
    lengths = paired(kernels, 2 * den, ref, 1, errors["reference"], effect_error)
    pairs = list(combinations(range(len(ports)), 2))
    # Six nearly equal geodesic distances in a round 3-sphere yield chords in
    # four Euclidean dimensions with an equal or larger shortest/longest ratio.
    measured = [
        lengths[pairs.index(pair)] for pair in combinations((0, 3, 4, 5, 6, 7), 2)
    ]
    low = min(x.a for x in measured)
    high = max(x.b for x in measured)
    quality = (low - iv.sqrt(iv.mpf(2) / 3) * high) / (1 + iv.sqrt(iv.mpf(2) / 3))
    assert quality.a > 0.30
    select = [pairs.index(pair) for pair in combinations(range(5), 2)]
    kernels = kernels[select]
    lengths = [lengths[i] for i in select]
    variations = {
        k: paired(kernels, 2 * den, v, 0, derr[k], effect_error)
        for k, v in changes.items()
    }
    dm = iv.matrix(5)
    for value, (i, j) in zip(lengths, combinations(range(5), 2), strict=True):
        dm[i, j] = dm[j, i] = value
    shifted = iv.matrix(4)
    for i in range(4):
        for j in range(4):
            shifted[i, j] = (dm[4, i] ** 2 + dm[4, j] ** 2 - dm[i, j] ** 2) / 2 - 3 * (
                i == j
            )
    for n in range(1, 5):
        assert iv.det(shifted[:n, :n]).a > 0
    a, v, w = geometry_jets(lengths, variations)
    ch = {
        "incoming": {
            "ordinary": iv.ln(2) / 3,
            "hidden": iv.mpf(source_bounds["hidden_modular_nats"]),
            "decorrelated": iv.mpf(source_bounds["decorrelation_modular_nats"]),
            "left": iv.mpf(0),
            "right": iv.mpf(0),
        },
        "retained": {
            k: iv.mpf(bounds) for k, bounds in material["modular_slopes"].items()
        },
    }
    methods = {}
    for adjuster in ("left", "right", "radius"):
        if adjuster == "radius":
            alpha = [iv.mpf(0)] * 3
            area = [g / 3 for g in w.g[:3]]
        else:
            i = KEYS.index(adjuster)
            assert v.g[i].a > 0 or v.g[i].b < 0
            alpha = [-g / v.g[i] for g in v.g[:3]]
            area = [g + angle * a.g[i] for g, angle in zip(a.g[:3], alpha, strict=True)]
        record = {
            "signed_area": {
                k: endpoints(x) for k, x in zip(KEYS[:3], area, strict=True)
            }
        }
        for name, charge in ch.items():
            kval = [
                charge[k] + (angle * charge[adjuster] if adjuster != "radius" else 0)
                for k, angle in zip(KEYS[:3], alpha, strict=True)
            ]
            eta = -kval[0] / area[0]
            record[name] = dict(
                ordinary_eta=endpoints(eta),
                residuals={
                    k: endpoints(z + eta * y)
                    for k, z, y in zip(KEYS[1:3], kval[1:], area[1:], strict=True)
                },
                predicted_area={
                    k: endpoints(-z / eta)
                    for k, z in zip(KEYS[1:3], kval[1:], strict=True)
                },
            )
        methods[adjuster] = record
    det = a.g[3] * v.g[4] - a.g[4] * v.g[3]
    assert det.b < 0
    simultaneous = {}
    for j, key in enumerate(KEYS[:3]):
        left = (-a.g[j] * v.g[4] + a.g[4] * v.g[j]) / det
        right = (-a.g[3] * v.g[j] + a.g[j] * v.g[3]) / det
        fixed_angle = -w.g[j] / w.g[4]
        simultaneous[key] = dict(
            fixed_outer_area_and_volume=dict(
                rotation_left=endpoints(left),
                rotation_right=endpoints(right),
                modular={
                    name: endpoints(c[key] + left * c["left"] + right * c["right"])
                    for name, c in ch.items()
                },
            ),
            fixed_cone_angle=dict(
                rotation_right=endpoints(fixed_angle),
                modular={
                    name: endpoints(c[key] + fixed_angle * c["right"])
                    for name, c in ch.items()
                },
            ),
        )
    return dict(
        stations=list(ports[:5]),
        apex_gram_eigenvalue_lower=3,
        all_radius_distance_error_lower=endpoints(quality)[0],
        pooled_range_effect_error=float(effect_error),
        distances=[endpoints(d) for d in lengths],
        outer_area=endpoints(a.v),
        interior_volume=endpoints(v.v),
        solid_angle=endpoints(w.v),
        methods=methods,
        two_rotation_jacobian=endpoints(det),
        fixed_geometry=simultaneous,
        charge_determinant_at_fixed_geometry=charge_rank_at_matched_geometry(
            simultaneous
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(certificate(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (live_survey.HERE / "live_bounds.json").write_text(data)
    else:
        print(data, end="")
