"""Four-station volume matching on the longer, fourth-order threshold relay graph."""


import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from .. import live_source, live_survey
from ..source import measure as source
from ....source_implication.common_front import measure as front
from ....source_implication.geometry_gain import simplex_area_volume


def read(board, stations, op, index, edges, state):
    c = np.einsum("kij,ji->k", op, state).real
    d, _, _ = front.graph(board.n, stations, (4 * c ** (-1 / 18))[index], edges)
    return np.asarray(simplex_area_volume(d))


def match(board, stations, op, index, edges, ref, changes):
    stations = tuple(map(int, stations[:4]))
    rho, product = source.form(board, source.formation.region(board), -1)
    g = live_survey.rotations(rho)["right"]
    channel = live_source.head_channel()
    val, vec = np.linalg.eigh(ref)
    k = -(vec * np.log(val)) @ vec.conj().T
    base = read(board, stations, op, index, edges, ref)

    def measure(state):
        return read(board, stations, op, index, edges, state)

    # Recompute every least path for every preparation, including finite-angle rotations.
    h = 1e-4
    derivatives = {
        tag: (measure(ref + h * changes[tag]) - measure(ref - h * changes[tag]))
        / (2 * h)
        for tag in ("ordinary", "decorrelated", "right")
    }
    rows = []
    for t in (1e-4, 1e-5, 1e-6):
        matches = {}
        for tag in ("ordinary", "decorrelated"):
            theta_slope = -derivatives[tag][1] / derivatives["right"][1]
            weights = (1 - t) * np.array(source.RAILS["reference"]) + t * np.array(
                source.RAILS[tag]
            )
            heads = (1 - t) * rho + t * product if tag == "decorrelated" else rho

            def prepare(theta):
                u = expm(-1j * theta * g)
                return live_source.apply_channel(
                    channel, u @ heads @ u.conj().T, weights
                )

            def residual(theta):
                return measure(prepare(theta))[1] - base[1]

            a, b = sorted([0.0, float(2 * t * theta_slope)])
            if residual(a) * residual(b) > 0:
                matches[tag] = dict(bracketed=False)
                continue
            theta = brentq(residual, a, b, xtol=1e-13)
            state = prepare(theta)
            a1, v1 = measure(state)
            vals = np.linalg.eigvalsh(state)
            matches[tag] = dict(
                bracketed=True,
                rotation=theta,
                volume_error=float(v1 - base[1]),
                area_secant=float((a1 - base[0]) / t),
                modular_secant=float(np.trace(k @ (state - ref)).real / t),
                entropy_secant=float(
                    (-np.dot(vals, np.log(vals)) + np.dot(val, np.log(val))) / t
                ),
            )
        row = dict(fraction=t, arms=matches)
        if all(m["bracketed"] for m in matches.values()):
            ordinary, held = matches["ordinary"], matches["decorrelated"]
            eta = -ordinary["modular_secant"] / ordinary["area_secant"]
            row.update(
                ordinary_gain=eta,
                held_residual=held["modular_secant"] + eta * held["area_secant"],
            )
        rows.append(row)
    return dict(
        stations=list(stations),
        reference_area=float(base[0]),
        reference_volume=float(base[1]),
        finite_steps=rows,
    )
