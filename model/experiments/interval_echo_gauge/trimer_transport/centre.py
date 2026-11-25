"""Read a constituent and its complementary colour phase on a single material clock."""

import argparse
import json
from functools import cache
from math import log, pi, sqrt

import networkx as nx
import numpy as np
from scipy.linalg import expm

from . import measure

ANGLES = (pi / 4, 2 * pi / 3, pi)
PICKS = ((0,), (1, 2), (0, 1, 2))
WEIGHTS = (1, 1, -2)


def endpoint_word(phase):
    """Two addressed private rails, including the vacant and doubled binary words."""
    K = measure.KAPPA
    t = (phase % (2 * pi)) / (2 * K + 1)
    x, c = measure.native_law.local_term_parts()
    labs = [measure.native_law.SRC.index(w) for w in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    h = (x + K * c)[np.ix_(labs, labs)]
    na = np.array([0.0, 0.0, 1.0, 1.0])
    nb = np.array([0.0, 1.0, 0.0, 1.0])
    z = np.diag(np.exp(-1j * K * (pi / K) * na))
    half = expm(-1j * t * h / 2)
    correction = (-t * (K + 0.5)) % (2 * pi)
    unitary = np.diag(np.exp(-1j * correction * (na + nb))) @ z @ half @ z @ half
    np.testing.assert_allclose(
        unitary, np.diag([1, 1, 1, np.exp(1j * phase)]), rtol=0, atol=3e-14
    )
    return np.diag(unitary)


def routing(routers):
    board = measure.carrier.federation(0)
    graph = nx.Graph()
    graph.add_edges_from(board.seams)
    # Two initialized spare home ports contain a single clock load.
    clock = next(
        (p, q)
        for p, q in board.seams
        if p // 12 == q // 12 == 0
        and p not in sum((s["source"] for s in routers), [])
        and q not in sum((s["source"] for s in routers), [])
        and board.degree[p] == board.degree[q] == 5
    )
    # After routing the source, the original home is empty on the admitted code.
    # To address a remote rail, transpose only the active clock address with a
    # vacant, private neighbour of that rail and invert the same permutation.
    result = []
    for receiver in routers:
        for target in receiver["triangle"]:
            slots = sorted(
                s
                for s in board.neighbours[target]
                if board.degree[s] == 5 and s not in receiver["triangle"]
            )
            neighbour = int(slots[0])
            path = min(nx.all_shortest_paths(graph, clock[0], neighbour))
            walk = list(zip(path, path[1:]))
            result.append(
                dict(
                    target=target,
                    ancilla=neighbour,
                    swaps=[list(e) for e in walk + walk[:-1][::-1]],
                )
            )
    return list(clock), result


@cache
def clock_instruments():
    contacts, private = measure.operators(measure.KAPPA)
    h = contacts[0, 1][:2, :2]
    phase = private[1][:2, :2]
    forward = expm(-1j * (pi / 2) * h)
    backward = expm(-1j * (3 * pi / 2) * h)
    minus_i = expm(-1j * (pi / (2 * measure.KAPPA)) * phase)
    plus_i = expm(-1j * (3 * pi / (2 * measure.KAPPA)) * phase)
    return minus_i @ forward[:, 0], (backward @ plus_i, backward)


def read_clock(density):
    # An inverse half contact and (for X) a private phase convert quadratures
    # to final ordinary position counts, without an abstract X/Y detector.
    coordinates = []
    for reader in clock_instruments()[1]:
        sharp = (reader @ density @ reader.conj().T).diagonal().real
        coordinates.append(float(sharp[0] - sharp[1]))
    return complex(*coordinates)


def fringe(tensor, theta, chosen, weights):
    prob = abs(tensor) ** 2
    prob = prob / prob.sum()
    clock, _ = clock_instruments()
    result = np.zeros((2, 2), complex)
    words = [endpoint_word(theta * w) for w in weights]
    for ind in measure.TUPLES:
        arms = clock.copy()
        for i in chosen:
            for rail, gate in enumerate(words):
                occupied = int(ind[i] == rail)
                arms *= gate[[occupied, 2 + occupied]]
        result += prob[ind] * np.outer(arms, arms.conj())
    return read_clock(result)


def model(theta, pick):
    single = (2 * np.exp(1j * theta) + np.exp(-2j * theta)) / 3
    return {1: single, 2: single.conjugate(), 3: 1 + 0j}[len(pick)]


def run():
    sec, arms, _, _, _, _, finite, src = measure.trimer.source()
    rails = np.delete(arms, 3, axis=1)
    clock, routes = routing(measure.geometry(rails))
    state = measure.tensor(sec, rails, finite)
    admit = float(np.vdot(state, state).real)
    source = src["filter_acceptances"][-1]
    rng = np.random.default_rng(12113571)
    rows = []
    for theta in ANGLES:
        for pick in PICKS:
            for name, weights in (
                ("cartan", WEIGHTS),
                ("number", (1, 1, 1)),
                ("parked", (0, 0, 0)),
            ):
                z = fringe(state, theta, pick, weights)
                record = []
                for component in (z.real, z.imag):
                    bright = (1 + component) / 2
                    record.append(
                        measure.count(
                            rng,
                            measure.STARTS,
                            [
                                1 - source,
                                source * (1 - admit),
                                source * admit * bright,
                                source * admit * (1 - bright),
                            ],
                        )
                    )
                readings = [(r[2] - r[3]) / (r[2] + r[3]) for r in record]
                radii = [sqrt(2 * log(2 * 128 / 0.01) / (r[2] + r[3])) for r in record]
                target = model(theta, pick)
                if name == "number":
                    target = np.exp(1j * theta * len(pick))
                if name == "parked":
                    target = 1 + 0j
                fits = all(
                    abs(value - forecast) <= width
                    for value, forecast, width in zip(
                        readings, (target.real, target.imag), radii, strict=True
                    )
                )
                rows.append(
                    dict(
                        theta=theta,
                        chosen=list(pick),
                        word=name,
                        counts=record,
                        readings=readings,
                        radii=radii,
                        fits=fits,
                    )
                )
    standards = []
    for theta in ANGLES:
        for colour in range(3):
            # One sharp material load on the selected address; all steps use
            # the same actual controlled native column as the trimer.
            gate = endpoint_word(theta * WEIGHTS[colour])
            initial, _ = clock_instruments()
            arms = initial * gate[[1, 3]]
            z = read_clock(np.outer(arms, arms.conj()))
            both = []
            for component in (z.real, z.imag):
                both.append(
                    measure.count(
                        rng, measure.CAL, [(1 + component) / 2, (1 - component) / 2]
                    )
                )
            standards.append(dict(theta=theta, colour=colour, counts=both))
    return dict(
        clock=clock,
        routes=routes,
        starts=measure.STARTS,
        rows=rows,
        single_load_standards=standards,
        calibration_starts=measure.CAL,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    output = measure.HERE / "centre.json"
    if args.write:
        output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(data, json.loads(output.read_bytes()), "centre")
    print(len(data["rows"]), "clock settings")
