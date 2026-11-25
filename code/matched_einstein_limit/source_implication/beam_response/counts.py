"""Independent synthetic five-pixel counts, including the outside outcome."""

import argparse
import json
import numpy as np
from mpmath import iv

from . import certificate, measure


def draw(rng, n, probabilities):
    # Each fixed-size batch and each entry of it fit into a signed machine integer.
    chunk = min(n, 2**58)
    assert n % chunk == 0
    rows = []
    for hits in probabilities.T:
        probs = np.r_[hits, 1 - hits.sum()]
        assert np.min(probs) >= 0
        totals = (
            rng.multinomial(chunk, probs, size=n // chunk).astype(object).sum(axis=0)
        )
        rows.append([int(c) for c in totals])
    return rows


def area_intervals(counts, n):
    # Each category in one multinomial is marginally binomial. A union bound
    # covers both species, all five preparations, five inputs and five hit pixels.
    # An additional factor two covers both tails. The miss counts are also archived.
    s = iv.ln(2 * 2 * 5 * 5 * 5 / iv.mpf(".01"))
    x = (
        np.array(
            [[iv.mpf(c) for c in row[:5]] for row in counts["equilibrium"]],
            dtype=object,
        ).T
        / n
    )

    def nominal(z):
        return float(z.mid)

    center = np.vectorize(nominal)(x)

    def axes(w):
        return np.c_[np.ones(5), measure.BEARINGS - w @ measure.BEARINGS / center.sum()]

    out_axis, inp_axis = axes(center.sum(1)), axes(center.sum(0))

    # These fixed translations preserve the determinant equation. Convert actual
    # integer-bearing / UNIT coordinates and binary-floating translations exactly.
    def exact_axis(a):
        weights = np.array(
            [
                [iv.mpf(1), *(iv.mpf(int(k)) / measure.UNIT for k in row)]
                for row in measure.PIXELS
            ],
            dtype=object,
        )
        shifts = [iv.mpf(0), iv.mpf(float(1 - a[0, 1])), iv.mpf(float(-a[0, 2]))]
        return weights - np.array(shifts, dtype=object)

    output = exact_axis(out_axis)
    inputs = exact_axis(inp_axis)

    def read(rows):
        p = np.empty((5, 5), dtype=object)
        for i, row in enumerate(rows):
            assert sum(row) == n
            for j, num in enumerate(row[:5]):
                radius = iv.sqrt(2 * iv.mpf(num) * s) + 2 * s
                p[j, i] = (iv.mpf(num) + iv.mpf([-1, 1]) * radius) / n
        return output.T @ p @ inputs

    moments = {k: read(data) for k, data in counts.items()}
    ref = moments["equilibrium"]
    co = np.array(
        [[certificate.cofactor(ref, i, j) for j in range(3)] for i in range(3)],
        dtype=object,
    )
    det = sum(a * b for a, b in zip(co[0], ref[0]))
    p = ref[0, 0]
    co[0, 0] -= 3 * det / p

    def slope(d):
        return 4 * sum((co * d).ravel()) / p**3

    rotation = slope((moments["plus"] - moments["minus"]) / 2)
    alpha = json.loads(
        (measure.HERE.parent / "null_optics/certificate.json").read_text()
    )["matched_pairing"]
    answer = {}
    for name, k in [
        ("ordinary", "ordinary_angle_interval"),
        ("pairing", "angle_interval"),
    ]:
        a, b = alpha[k]
        angle = iv.mpf([int(np.floor(a)), int(np.ceil(b))])
        answer[name] = slope(moments[name] - ref) + angle * rotation
    return {k: str(v) for k, v in answer.items()}


def detector_laws():
    f = measure.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = measure.screens(board, route)
    states = measure.aperture.survey.preparations()
    laws = {}
    for q in (-1, 1):
        table = measure.effects(board, route, rings, q, 4.0)
        laws[str(q)] = {
            name: measure.probabilities(table, states[name])
            for name in ("equilibrium", "ordinary", "pairing", "plus", "minus")
        }
    return laws


def run(n=2**72, laws=None):
    iv.dps = 30
    rng = np.random.default_rng(638512)
    data = {}
    for q, settings in (detector_laws() if laws is None else laws).items():
        counts = {
            name: draw(rng, n, probabilities)
            for name, probabilities in settings.items()
        }
        data[str(q)] = dict(
            hit_pixels_then_miss=counts, signed_slopes_99pct=area_intervals(counts, n)
        )
    return dict(dwell=4, shots_per_pixel_per_preparation=n, seed=638512, species=data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if opt.write:
        (measure.HERE / "counts.json").write_text(result)
    else:
        print(result, end="")
