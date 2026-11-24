"""Measure the one-anchor consequences of a compact vacancy on the pinned law."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import (
    defect_generator,
    free_seam,
)
from model.experiments.material_inertia import neutral_defect
from . import experiment
from .gateway_acquisition import prepare_layout

K = candidate.KAPPA
SHOTS = 2**20


def preparation(board):
    charts = experiment.gateway_table(board)
    local, antipodes = experiment.geometry()
    d = shortest_path(
        board.adjacency(), directed=False, indices=np.arange(12), unweighted=True
    )
    same = [i for i, chart in enumerate(charts) if chart == charts[0]]
    distance = d.min(axis=0).reshape(-1, 12).min(axis=1)
    away = max(same, key=lambda i: (distance[i], -i))
    tips = charts[0]
    avoid = {p for pair in antipodes if set(pair).intersection(tips) for p in pair}
    launch = min(set(range(local.n)) - avoid)
    prep = prepare_layout(tips)
    accepted = prep.accepted_vectors[:, launch]
    success = np.linalg.norm(accepted) ** 2
    fails = prep.first_click[:, :, launch].sum(axis=1)
    np.testing.assert_allclose(fails.sum() + success, 1, atol=1e-11)
    # Both signs see the same isolated pulse on the unit-filled medium.
    law = candidate.native_law
    seam = free_seam(K)
    for sign in (-1, 1):
        labels = ((1, 1), (1 + sign, 1), (1, 1 + sign))
        indices = [law.SRC.index(pair) for pair in labels]
        np.testing.assert_allclose(
            seam[np.ix_(indices, indices)],
            [[0, 0, 0], [0, K + 0.5, -0.5], [0, -0.5, K + 0.5]],
        )
        blank = np.eye(law.DIM)[:, indices]
        echo = law.echo()
        for outcome in (-1, 0, 1):
            tags = np.repeat([False, *(v == outcome for _, v in law.LABELS)], law.N_SRC)
            np.testing.assert_allclose(
                echo.T @ (tags[:, None] * (echo @ blank)),
                blank * [a - b == outcome for a, b in labels],
                atol=3e-14,
            )
    states = []
    for i in (0, away):
        v = np.zeros(board.n, complex)
        v[12 * i : 12 * i + 12] = accepted / np.sqrt(success)
        states.append(v)
    return (
        states,
        success,
        {
            "cells": [0, away],
            "minimum_port_distance": int(distance[away]),
            "launch_local_port": launch,
            "single_attempt_acceptance": float(success),
            "first_gateway_failure_counts_then_accepted": np.random.default_rng(83019)
            .multinomial(SHOTS, [*fails, success])
            .tolist(),
        },
    )


def read_energy(h, v, success, rng):
    # Uniform host energy is subtracted during assembly. Unequal-load arms require a reservoir.
    action = h @ v
    energy = np.vdot(v, action).real
    residual = np.linalg.norm(action - energy * v)
    if residual > 1e-10:
        raise ArithmeticError(
            "The monitored state does not support a single stationary line"
        )
    t0 = 1 / 128
    test = expm_multiply((-1j * t0) * h, v)
    np.testing.assert_allclose(np.vdot(v, test), np.exp(-1j * t0 * energy), atol=1e-12)
    # Larger waits use the certified one-column reduction, with error <= time * residual.
    phase = 0.0
    rows = []
    for k in range(18):
        t = t0 * 2**k
        xy = np.array([np.cos(energy * t), np.sin(energy * t)])
        xy = np.clip(np.round(xy, 12), -1, 1)
        counts = [
            rng.multinomial(
                SHOTS, [success * (1 + x) / 2, success * (1 - x) / 2, 1 - success]
            )
            for x in xy
        ]
        estimate = complex(*[(a - b) / (a + b) for a, b, _ in counts])
        arg = np.angle(estimate)
        phase = arg + 2 * np.pi * round((2 * phase - arg) / (2 * np.pi))
        rows.append([t, [r.tolist() for r in counts]])
    return {
        "energy_above_filled_host": float(energy),
        "eigen_residual": float(residual),
        "ramsey_estimate": float(phase / t),
        "counts_by_wait": rows,
    }


def run():
    board = federation(0)
    (left, right), success, prep = preparation(board)
    rng = np.random.default_rng(83020)
    pairs = PairSector(board)
    a, b = pairs.aa, pairs.bb
    vv = left[a] * right[b] + right[a] * left[b]
    x, y, adj, exchange, native_shifted, _ = neutral_defect.generators(board)
    ev = right[x] * left[y]
    for v in (left, right, vv, ev):
        np.testing.assert_allclose(np.linalg.norm(v), 1, atol=1e-13)
    outcomes = {}
    for ref in (0, 1):
        vacancy, excess = [defect_generator(board, sign, K, ref) for sign in (-1, 1)]
        vacancy_cost = K * (6 - 2 * (1 - ref) * (6 - board.degree))
        double = pairs.exchange() + sparse.diags(
            vacancy_cost[a] + vacancy_cost[b] - 2 * K * pairs.adjacent
        )
        neutral = exchange + sparse.diags(
            12 * K + 2 * K * ((1 - ref) * (board.degree[y] - board.degree[x]) + adj)
        )
        if ref == 0:
            diff = neutral - native_shifted - 12 * K * sparse.eye(len(ev))
            np.testing.assert_allclose(diff.data, 0, atol=2e-13)
        first, second = board.seams[0]
        adjacent_excess = (
            double[pairs.index(first, second), pairs.index(first, second)]
            - vacancy[first, first]
            - vacancy[second, second]
        )
        outcomes[str(ref)] = {
            "one_vacancy": read_energy(vacancy, left, success, rng),
            "one_excess": read_energy(excess, right, success, rng),
            "two_vacancies": read_energy(double, vv, success**2, rng),
            "neutral_pair": read_energy(neutral, ev, success**2, rng),
            "adjacent_sharp_two_vacancies_minus_singles": float(adjacent_excess),
        }
    native = outcomes["0"]
    light = native["one_vacancy"]["ramsey_estimate"]
    energies = {
        k: float(511 * v["ramsey_estimate"] / light)
        for k, v in native.items()
        if isinstance(v, dict)
    }
    bounds = []
    for level in (0, 1, 2):
        car = federation(level)
        h = defect_generator(car, -1, K, 0)
        radius = abs(h).sum(axis=1).A1 - abs(h.diagonal())
        bounds.append(
            {
                "cells": car.cells,
                "formation_energy_bounds": [
                    float(min(h.diagonal() - radius)),
                    float(max(h.diagonal() + radius)),
                ],
            }
        )
    return {
        "source": prep,
        "shots_per_quadrature": SHOTS,
        "seed": 83020,
        "records_by_degree_reference": outcomes,
        "native_anchored_keV_from_counts": energies,
        "all_single_vacancy_spectral_enclosures": bounds,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output is None:
        from ..electron_field.replay import compare

        stored = Path(__file__).with_suffix(".json")
        compare(result, json.loads(stored.read_text()), "anchored_vacancy")
        print(
            "Reproduced the anchored vacancy, both pair sectors and the reference control."
        )
    else:
        with args.output.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
