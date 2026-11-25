"""Independent finite-clock angular and dimer-count replay including every miss."""

import argparse
import json

import numpy as np
from mpmath import iv
from scipy.linalg import expm

from ..beam_response.counts import draw
from ..material_clocks.measure import blank_clock
from . import measure

HERE = measure.HERE


def area_band(rows, n, tail):
    # Independent multinomial rows. Bound all hit probabilities with two-sided Bernstein.
    p = np.empty((5, 5), dtype=object)
    for i, row in enumerate(rows):
        assert sum(row) == n
        for j, hits in enumerate(row[:5]):
            radius = iv.sqrt(2 * iv.mpf(hits) * tail) + 2 * tail
            p[j, i] = (iv.mpf(hits) + iv.mpf([-1, 1]) * radius) / n
    bearings = np.array(
        [
            [iv.mpf(int(i)) / measure.beam.UNIT, iv.mpf(int(j)) / measure.beam.UNIT]
            for i, j in measure.beam.PIXELS
        ]
    )
    # Translating the bearings before centering reduces interval inflation.
    mid = np.array([[float(x.mid) for x in row] for row in p])

    def basis(axis):
        center = (mid.sum(axis=axis) @ measure.beam.BEARINGS) / mid.sum()
        shifted = bearings - np.array([iv.mpf(float(v)) for v in center])
        return np.c_[np.ones(5, dtype=int), shifted]

    m = basis(1).T @ p @ basis(0)
    c = m[1:, 1:] - np.outer(m[1:, 0], m[0, 1:]) / m[0, 0]
    return 4 * (c[0, 0] * c[1, 1] - c[0, 1] * c[1, 0]) / m[0, 0] ** 2


def run(n=2**62):
    iv.dps = 32
    f = measure.beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = measure.beam.screens(board, route)
    x = measure.beam.aperture.survey.preparations()["equilibrium"]
    rng = np.random.default_rng(628531)
    tail = iv.ln(iv.mpf(2 * 2 * 5 * 5 * 3) / iv.mpf(".01"))
    results = {}
    for t in measure.TIMES:
        dimer = blank_clock(1, 8)
        # Compensate the known scalar energy. Detector reads the opposite dimer endpoint.
        pd = float(abs(expm(-1j * t * dimer)[1, 0]) ** 2)
        clockshots = 2**40
        clock = int(rng.binomial(clockshots, pd))
        result = {
            "clock_counts": [clock, clockshots - clock],
            "clock_probability": pd,
            "dwell_from_calibrated_clock": 2 * np.pi
            - 2 * np.arcsin(np.sqrt(clock / clockshots)),
        }
        for q in (-1, 1):
            p = measure.beam.probabilities(
                measure.beam.effects(board, route, rings, q, t), x
            )
            rows = draw(rng, n, p)
            a = area_band(rows, n, tail)
            assert a.b < 0 if (q == -1 and t == 4.5) else a.a > 0
            result[str(q)] = dict(counts_by_sender=rows, signed_area_99pct=str(a))
        results[str(t)] = result
    # Separate population replay samples a centered stencil of real clock appointments.
    stencils = []
    for q in (-1, 1):

        def read(t):
            return measure.beam.image(
                measure.beam.probabilities(
                    measure.beam.effects(board, route, rings, q, t), x
                )
            )[0]

        base = read(4)
        s = []
        for h in (1 / 32, 1 / 64, 1 / 128):
            plus, minus = read(4 + h), read(4 - h)
            _, b, c = measure.optical_matrices(
                [base, (plus - minus) / (2 * h), (plus + minus - 2 * base) / h**2]
            )
            s.append(
                dict(
                    half_tick=h,
                    expansion=float(np.trace(b)),
                    clock_focusing=-float(np.trace(c)),
                )
            )
        control = measure.signal_jets(board, route, rings, q, 4, strength=0)
        bmat = np.c_[np.ones(5), measure.beam.BEARINGS]
        probs = np.array([measure.beam.probabilities(v, x) for v in control])
        det, b, c = measure.optical_matrices(measure.matrix_jets(bmat.T @ probs @ bmat))
        stencils.append(
            dict(
                charge=q,
                stencils=s,
                zero_contrast_dwell_four=dict(
                    area=det,
                    expansion=float(b.trace()),
                    clock_focusing=-float(c.trace()),
                ),
            )
        )
    return dict(
        seed=628531,
        shots_per_pixel=n,
        clock_branch="pi < t < 2 pi",
        records=results,
        population_comparisons=stencils,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "counts.json").write_text(output)
    else:
        print(output, end="")
