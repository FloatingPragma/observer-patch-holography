"""Measure single-arm and two-arm transport using the energy-selected native trimer."""

import argparse
import json
from itertools import combinations, product
from math import log, pi, sqrt
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.candidate_laws.interval_echo.openwave.row_13_baryon.arm_symmetry import (
    measure as trimer,
)
from model.executable_model import carrier

HERE = Path(__file__).resolve().parent
STARTS = 2**32
CAL = 2**18
MENUS = (("AB", "AB"), ("AB", "BA"), ("BAB", "AB"), ("A", "B"))
TUPLES = list(product(range(3), repeat=3))


def geometry(arms):
    board = carrier.federation(0)
    graph = nx.Graph()
    graph.add_edges_from(board.seams)
    owners = sorted(
        {p // 12 for e in board.external if 0 in {q // 12 for q in e} for p in e} - {0}
    )
    routes = []
    for original, owner in zip(arms, owners, strict=True):
        region = tuple(range(12 * owner, 12 * owner + 12))
        triangle = next(
            t
            for t in combinations(region, 3)
            if all(board.degree[j] == 5 for j in t)
            and all(graph.has_edge(i, j) for i, j in combinations(t, 2))
        )
        local = graph.subgraph((*range(12), *region))
        moves = []
        for source, destination in zip(original, triangle, strict=True):
            path = min(nx.all_shortest_paths(local, source, destination))
            forward = list(zip(path, path[1:]))
            moves.extend(forward + forward[:-1][::-1])
        routes.append(
            dict(
                owner=owner,
                source=list(map(int, original)),
                triangle=list(triangle),
                pulses=[list(map(int, e)) for e in moves],
            )
        )
    return routes


def operators(kappa):
    ex, co = native_law.local_term_parts()
    h = (ex + kappa * co)[:9, :9]
    loads = np.eye(3, dtype=int)
    onsite = [np.diag(kappa * loads[:, p] ** 2) for p in range(3)]
    contacts = {}
    for a, b in combinations(range(3), 2):
        idx = [native_law.SRC.index((w[a], w[b])) for w in loads]
        contacts[a, b] = h[np.ix_(idx, idx)]
    return contacts, onsite


def native_word(word, kappa=KAPPA, operation="forward"):
    contacts, private = operators(kappa)
    pulses = []
    for letter in word:
        a, b = (0, 1) if letter == "A" else (1, 2)
        t = pi / 2
        # Cancel the common κ+1/2 phase on the two single-load addresses.
        # This exposure uses the original clock calibration on the control as well.
        correction = (-t * (KAPPA + 0.5)) % (2 * pi)
        pulses.extend(
            [
                (contacts[a, b], t),
                (private[a], correction / KAPPA),
                (private[b], correction / KAPPA),
            ]
        )
        if letter == "B":
            pulses.append((private[2], (pi / 2) / KAPPA))
    if operation in ("inverse", "transpose"):
        pulses.reverse()
    answer = np.eye(3, dtype=complex)
    for generator, t in pulses:
        # Native blank contact and private spectra are integral for both laws.
        wait = (-t) % (2 * pi) if operation == "inverse" else t
        answer = expm(-1j * wait * generator) @ answer
    return answer


def tensor(sec, arms, vector):
    return vector[[sec.index[tuple(sorted(w))] for w in product(*arms)]].reshape(
        3, 3, 3
    )


def detector(state, first, left, right, dephase=False):
    if dephase:
        return np.einsum(
            "ai,bj,ck,ijk->abc",
            abs(first) ** 2,
            abs(left) ** 2,
            abs(right) ** 2,
            abs(state) ** 2,
        )
    return abs(np.einsum("ai,bj,ck,ijk->abc", first, left, right, state)) ** 2


def count(rng, starts, prob):
    p = np.maximum(0, np.round(prob, 14))
    if abs(p.sum() - 1) > 2e-12:
        raise ValueError("Missing probability")
    return rng.multinomial(starts, p / p.sum()).tolist()


def expected(cal):
    table = np.asarray(cal)
    if (
        table.shape != (3, 3)
        or min(table.ravel()) < 0
        or np.any(table.sum(axis=1) != CAL)
    ):
        raise ValueError("Invalid calibration")
    return np.array(
        [table[3 - b - c, a] / (6 * CAL) if b != c else 0.0 for a, b, c in TUPLES]
    )


def score(counts, cal):
    if len(counts) != 29 or min(counts) < 0 or sum(counts) != STARTS:
        raise ValueError("Incomplete source records")
    admitted = sum(counts[2:])
    if admitted == 0:
        return dict(admissions=0, covered=None)
    measured = np.array(counts[2:]) / admitted
    error = max(abs(measured - expected(cal)))
    # At most 4096 frequencies across all complete held and calibration banks.
    factor = log(2 * 4096 / 0.01) / 2
    allowance = sqrt(factor / admitted) + sqrt(factor / CAL) / 6
    return dict(
        admissions=admitted,
        max_discrepancy=float(error),
        allowance=allowance,
        covered=bool(error <= allowance),
    )


def run():
    sec, arms, _, _, _, _, finite, source = trimer.source()
    rails = np.delete(arms, 3, axis=1)
    routes = geometry(rails)
    source_rate = source["filter_acceptances"][-1]
    admitted = tensor(sec, rails, finite)
    _, _, _, _, gateway_count = trimer.apparatus()
    en, vec = np.linalg.eigh(sec.hamiltonian(0, False) + np.diag(gateway_count / 2))
    ctrl, rates = trimer.filters(en, vec, trimer.TARGET, sec.index[(1, 8, 9)])
    control = tensor(sec, rails, ctrl)
    rng = np.random.default_rng(12113548)
    banks = []
    for x, y in MENUS:
        first = native_word(x)
        other = native_word(y)
        cal = {}
        for direction in ("inverse", "transpose"):
            transfer = first @ native_word(y, operation=direction)
            cal[direction] = [count(rng, CAL, abs(col) ** 2) for col in transfer.T]
        settings = {}
        for case in ("native", "dephased", "one_route_parked", "zero_contrast"):
            if case == "zero_contrast":
                rate = rates[-1]
                signal = detector(
                    control, native_word(x, 0), native_word(y, 0), native_word(y, 0)
                )
            else:
                rate = source_rate
                right = np.eye(3) if case == "one_route_parked" else other
                signal = detector(admitted, first, other, right, case == "dephased")
            bins = count(
                rng,
                STARTS,
                np.r_[1 - rate, rate * (1 - signal.sum()), rate * signal.ravel()],
            )
            settings[case] = dict(
                counts=bins,
                scores={mode: score(bins, reads) for mode, reads in cal.items()},
            )
        banks.append(dict(first=x, pair=y, calibration=cal, records=settings))
    return dict(
        starts=STARTS,
        calibration_starts=CAL,
        reference_port_owner=0,
        filters=source["filter_acceptances"],
        finite_source_vector_error=source["finite_to_ideal_vector_upper"],
        zero_contrast_filter=rates,
        admitted_given_source=float(np.vdot(admitted, admitted).real),
        routes=routes,
        banks=banks,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    raw = json.dumps(data, indent=2, allow_nan=False) + "\n"
    f = HERE / "counts.json"
    if args.write:
        f.write_text(raw)
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(data, json.loads(f.read_bytes()), "trimer_transport")
    for row in data["banks"]:
        print(
            row["first"],
            row["pair"],
            {
                k: [v["scores"][name]["covered"] for name in ("inverse", "transpose")]
                for k, v in row["records"].items()
            },
        )
