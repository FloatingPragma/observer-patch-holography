"""Leave every sharp emission and free history fixed. Change only the final receiver optics."""

import argparse
import json
from itertools import combinations

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response.four_ticks import parse_interval
from . import certificate as cert
from . import measure


def odd_minors(settings):
    odd = {
        tag: [
            parse_interval(settings["forward"][tag][x])
            - parse_interval(settings["reverse"][tag][x])
            for x in ["expansion", "focusing"]
        ]
        for tag in measure.PHASES
    }
    minors = {}
    for tri in combinations(odd, 3):
        a, b, c = (odd[t] for t in tri)
        det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
        minors["/".join(tri)] = str(det)
    return minors


def receipt():
    iv.dps = 32
    f = measure.optics.beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = measure.optics.beam.screens(board, route)
    blanks = {tuple(sorted(ab)) for ab in board.seams} - {
        tuple(sorted(ab)) for ab in (route[:2], route[2:])
    }
    assert all(
        tuple(sorted(r[start : start + 2])) in blanks for r in rings for start in (0, 2)
    )
    state, rad = cert.sharp.bounded_sources(sp.eye(4) / 4)["state"]
    nominal = measure.optics.beam.aperture.survey.preparations()["equilibrium"]
    coeff = cert.sharp.prior.weights()[[7]]
    qs = {}
    for q in (-1, 1):
        waves = []
        exact = []
        for odd in (False, True):
            waves.append(measure.amplitudes(board, route, rings, q, odd, 8))
            h = cert.sharp.native.generator(board, route, q, odd)
            exact.append(
                np.array(cert.sharp.amplitude_jets(h, *rings, coeff), dtype=object)
            )
        ww = np.array(waves)
        settings = {}
        calibrated = {}
        for side in ["forward", "reverse"]:
            kernel = ww if side == "forward" else ww.transpose(0, 1, 3, 2, 5, 4)
            fixed = (
                exact
                if side == "forward"
                else [v.transpose(0, 1, 2, 4, 3, 6, 5) for v in exact]
            )
            tables = {}
            jac = {}
            for tag, (first, second) in measure.PHASES.items():
                # The two blank seams are disjoint. Relative phases between
                # their output blocks multiply terminal rows by harmless phases.
                g = measure.splitter(q, first)
                g[2:4, 2:4] = measure.splitter(q, second)[:2, :2]
                p = measure.probabilities(kernel, nominal, np.eye(5), g)
                transformed = []
                for v in fixed:
                    a, b = cert.split_axis(
                        v[:, 0], v[:, 1], 2, [(0, 1, first), (2, 3, second)]
                    )
                    jets = cert.sharp.probability_jets(np.stack((a, b), axis=1))[
                        :, :, 0
                    ]
                    transformed.append(jets)
                op = np.array(transformed, dtype=object) * sp.Rational(1, 25)
                theta, focus, p00, jac[tag] = cert.readings(op, state, rad)
                read = measure.read(p)
                assert (
                    read["expansion"] in theta
                    and read["focusing"] in focus
                    and p[0, 0, 0] in p00
                )
                # The receive basis spans the same subspace in every setting,
                # even conditional on the full random launch-pixel record.
                tables[tag] = dict(
                    expansion=str(theta),
                    focusing=str(focus),
                    pixel00=str(p00),
                    pixel22_nominal=float(p[0, 2, 2]),
                    hit_by_launch=p[0].sum(axis=0).tolist(),
                )
            for t in measure.PHASES:
                np.testing.assert_allclose(
                    tables[t]["hit_by_launch"],
                    tables["plus_plus"]["hit_by_launch"],
                    atol=1e-15,
                )
            settings[side] = tables
            calibrated[side] = cert.moving_screen_test(jac)
        qs[str(q)] = dict(
            settings=settings,
            odd_minors=odd_minors(settings),
            moving_screen_tests=calibrated,
        )
    return dict(dwell=4, detector_edges=[[r[:2], r[2:4]] for r in rings], charges=qs)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    out = json.dumps(receipt(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (measure.HERE / "delayed.json").write_text(out)
    else:
        print(out, end="")
