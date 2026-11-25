"""Independent ray tests on the same density for reciprocal coherent phase settings."""

import argparse
import json
from itertools import combinations

import numpy as np
import sympy as sp
from mpmath import iv

from .. import certificate as sharp
from .. import measure as optics
from . import measure

HERE = measure.HERE


def split_axis(real, imag, axis, pairs):
    a, b = np.moveaxis(real, axis, 0), np.moveaxis(imag, axis, 0)
    c, d = 5 * a, 5 * b
    for first, second, sign in pairs:
        for i, j in [(first, second), (second, first)]:
            c[i] = 4 * a[i] - sign * 3 * b[j]
            d[i] = 4 * b[i] + sign * 3 * a[j]
    return np.moveaxis(c, 0, axis), np.moveaxis(d, 0, axis)


def split(amps, left, right):
    # Each gate has denominator five, with axes (order,time,out,in,coin,coin).
    real, imag = split_axis(amps[:, 0], amps[:, 1], 2, [(0, 1, right)])
    real, imag = split_axis(real, imag, 3, [(0, 1, left)])
    return np.stack((real, imag), axis=1)


def physical_tables(waves, left, right):
    return np.array(
        [sharp.probability_jets(split(w, left, right))[:, :, 0] for w in waves],
        dtype=object,
    ) * sp.Rational(1, 625)


def readings(ops, state, rad):
    moments = sharp.enclose_moments(ops, state, rad)
    # A gate has maximum absolute row and column sum 7/5. Applying both
    # sides costs <2 in amplitude uncertainty, including time derivatives.
    # Add a second allowance twice as large as that in enclose_moments.
    for order, m in enumerate(moments):
        moments[order] = m + iv.mpf([-1, 1]) * iv.mpf(2000 * 64**order) / 10**19
    jac = optics.matrix_jets(moments)
    _, b, c = optics.optical_matrices(jac)
    probs = []
    for re, im in ops[:, 0, :, 0, 0]:
        probs.append((sp.Matrix(re) + sp.I * sp.Matrix(im)) / sharp.prior.FRAC**2)
    center = sum(sharp.full.inner(p, x) for p, x in zip(probs, state))
    norm = max(sharp.prior.one_norm(p) for p in probs)
    p = sharp.full.ri_band(center, norm * rad + sp.Rational("3e-18"))
    return b.trace(), -c.trace(), p, jac


def moving_screen_test(optics_by_tag):
    # For a common linear optical equation in *any* transported output frame
    # the acceleration row lies in the span of the two position and velocity rows.
    data = {}
    for component in (0, 1):
        full = []
        for tag in ("plus_plus", "minus_minus", "plus_minus"):
            j, d, e = optics_by_tag[tag]
            for col in (0, 1):
                full.append(
                    [d[0, col], d[1, col], j[0, col], j[1, col], e[component, col]]
                )
        full = np.array(full, dtype=object)
        tested = []
        for rows in combinations(range(6), 5):
            m = full[list(rows)]
            mid = np.array([[float(x.mid) for x in row] for row in m])
            try:
                inv = np.linalg.inv(mid)
            except np.linalg.LinAlgError:
                continue
            # Only the outward dyadic computation determines invertibility.
            scale = 2**32
            pre = np.array(
                [[iv.mpf(int(round(x * scale))) / scale for x in row] for row in inv]
            )
            residual = np.eye(5, dtype=int) - pre @ m
            upper = max(sum(abs(x) for x in row).b for row in residual)
            tested.append((upper, rows))
        upper, rows = min(tested, key=lambda attempt: float(attempt[0]))
        data[str(component + 1)] = dict(
            column_rows=list(rows),
            identity_error_infinity_upper=str(upper),
            certified_full_rank=bool(upper < 1),
        )
    return data


def run():
    iv.dps = 32
    f = optics.beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = optics.beam.screens(board, route)
    coeff = sharp.prior.weights()[[7]]
    state, rad = sharp.bounded_sources(sp.eye(4) / 4)["state"]
    results = {}
    for charge in (-1, 1):
        waves = []
        for odd in (False, True):
            h = sharp.native.generator(board, route, charge, odd)
            waves.append(np.array(sharp.amplitude_jets(h, *rings, coeff), dtype=object))
        rev = [v.transpose(0, 1, 2, 4, 3, 6, 5) for v in waves]
        tags = {}
        raw = {}
        panels = {"forward": {}, "reverse": {}}
        for tag, (left, right) in measure.PHASES.items():
            a = readings(physical_tables(waves, left, right), state, rad)
            b = readings(physical_tables(rev, right, left), state, rad)
            raw[tag] = (a[0] - b[0], a[1] - b[1])
            panels["forward"][tag] = a[3]
            panels["reverse"][tag] = b[3]
            tags[tag] = dict(
                odd_expansion=str(raw[tag][0]),
                odd_clock_focusing=str(raw[tag][1]),
                forward_pixel00=str(a[2]),
                reverse_pixel00=str(b[2]),
            )
        # Each stationary tag has identical native density and hence the same
        # value of any tag-independent regional source, not only the bridge.
        # Eliminate both the source and any calibration using three readings.
        a = raw["plus_plus"]
        b = raw["minus_minus"]
        c = raw["plus_minus"]
        slopes = [(a[1] - d[1]) / (d[0] - a[0]) for d in [b, c]]
        det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
        assert 0 not in det
        assert slopes[0].b < slopes[1].a or slopes[1].b < slopes[0].a
        results[str(charge)] = dict(
            tags=tags,
            orientation_minor=str(det),
            required_nonaffinities=[str(v) for v in slopes],
            moving_screen_tests={
                side: moving_screen_test(rows) for side, rows in panels.items()
            },
        )
    return dict(
        splitter_triplet=["plus_plus", "minus_minus", "plus_minus"], charges=results
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "certificate.json").write_text(data)
    else:
        print(data, end="")
