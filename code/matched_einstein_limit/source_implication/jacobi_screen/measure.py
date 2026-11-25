"""Two screened derivatives of the complete angular image with respect to native dwell."""

import argparse
from pathlib import Path
import json
from math import comb

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ..beam_response import measure as beam

HERE = Path(__file__).resolve().parent
TIMES = (4.0, 4.5, 5.0)


def signal_jets(board, route, rings, charge, time, strength=8):
    """The twenty launch columns give the squared amplitudes and their two derivatives."""
    rows, cols = beam.indices(rings[1]), beam.indices(rings[0])
    signals = []
    for odd in (False, True):
        h = (
            beam.front.native.generator(board, route, charge, odd) / 2
            if strength == 8
            else beam.front.optics.ranging.generator(
                board, route, charge, strength, odd
            )
        )
        launches = np.zeros((h.shape[0], len(cols)))
        launches[cols] = np.eye(len(cols))
        h = -1j * h
        wave = expm_multiply(time * h, launches, traceA=time * h.diagonal().sum())
        amplitudes = []
        for order in range(3):
            amplitudes.append(wave[rows].reshape(5, 4, 5, 4).transpose(0, 2, 1, 3))
            if order < 2:
                wave = h @ wave
        f, ft, ftt = amplitudes

        def adj(a):
            return a.conj().swapaxes(-1, -2)

        mixed = adj(f) @ ft
        second = adj(f) @ ftt + adj(ft) @ ft
        signals.append([adj(f) @ f, mixed + adj(mixed), second + adj(second)])
    return np.array(signals).swapaxes(0, 1)


def matrix_jets(moments):
    """Differentiate centered detector moments, including throughput and both marginals."""
    m, dt, d2t = moments
    w = m[0, 0]
    normalized = m / w
    first = (dt - normalized * dt[0, 0]) / w
    second = (d2t - 2 * first * dt[0, 0] - normalized * d2t[0, 0]) / w

    def prod(a, b):
        return np.outer(a[1:, 0], b[0, 1:])

    return [
        2 * (normalized[1:, 1:] - prod(normalized, normalized)),
        2 * (first[1:, 1:] - prod(first, normalized) - prod(normalized, first)),
        2
        * (
            second[1:, 1:]
            - prod(second, normalized)
            - 2 * prod(first, first)
            - prod(normalized, second)
        ),
    ]


def response_jets(m, dm):
    """First source variations of J, J-dot and J-double-dot, at a fixed clock."""
    norm, change, out = [], [], []
    for order in range(3):
        term = m[order].copy()
        for k in range(order):
            term = term - comb(order, k) * norm[k] * m[order - k][0, 0]
        norm.append(term / m[0][0, 0])
        dt = dm[order].copy()
        for k in range(order + 1):
            dt = dt - comb(order, k) * norm[k] * dm[order - k][0, 0]
            if k < order:
                dt = dt - comb(order, k) * change[k] * m[order - k][0, 0]
        change.append(dt / m[0][0, 0])
        image = change[order][1:, 1:].copy()
        for k in range(order + 1):
            image = image - comb(order, k) * (
                np.outer(change[k][1:, 0], norm[order - k][0, 1:])
                + np.outer(norm[k][1:, 0], change[order - k][0, 1:])
            )
        out.append(2 * image)
    return out


def deformation_response(j, dj):
    det, b, c = optical_matrices(j)
    inverse = np.array([[j[0][1, 1], -j[0][0, 1]], [-j[0][1, 0], j[0][0, 0]]]) / det
    return (dj[1] - b @ dj[0]) @ inverse, (dj[2] - c @ dj[0]) @ inverse


def optical_matrices(jacobians):
    j, jt, jtt = jacobians
    inverse = np.array([[j[1, 1], -j[0, 1]], [-j[1, 0], j[0, 0]]])
    det = j[0, 0] * j[1, 1] - j[1, 0] * j[0, 1]
    return det, (jt @ inverse) / det, (jtt @ inverse) / det


def run():
    f = beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = beam.screens(board, route)
    x = beam.aperture.survey.preparations()["equilibrium"]
    basis = np.c_[np.ones(5), beam.BEARINGS]
    rows = []
    for time in TIMES:
        panels = {}
        for q in (-1, 1):
            signals = signal_jets(board, route, rings, q, time)
            probs = np.array([beam.probabilities(s, x) for s in signals])
            moments = basis.T @ probs @ basis
            j = matrix_jets(moments)
            det, b, c = optical_matrices(j)
            expansion = np.trace(b)
            shear = (b + b.T) / 2 - expansion * np.eye(2) / 2
            twist = (b[1, 0] - b[0, 1]) / 2
            panels[str(q)] = dict(
                jacobian=[v.tolist() for v in j],
                signed_area=det,
                throughput=moments[0, 0, 0] / 5,
                optical_deformation=b.tolist(),
                acceleration_map=c.tolist(),
                expansion=expansion,
                shear=shear.tolist(),
                twist=twist,
                trace_focusing_if_affine=-float(np.trace(c)),
            )
        rows.append(dict(dwell=time, charges=panels))
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    s = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(s)
    else:
        print(s, end="")
