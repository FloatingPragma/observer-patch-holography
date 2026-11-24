"""Two local contact histories and two enclosing cuts on compact native sources."""

import argparse
from functools import cache
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse import eye
import sympy as sp

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import CELL_SEAMS, federation
from model.experiments.interval_echo_encounter.conjugation import (
    defect_generator,
    free_seam,
)
from model.experiments.material_inertia.pair_spectrum import probe

HERE = Path(__file__).resolve().parent
SHOTS = 2**22
SEED = 492301
K = candidate.KAPPA
SCORE = (probe.BINS == 1).astype(int) - (probe.BINS == -1)
PRECISION_SHOTS = 2**28


def simultaneous_radius(shots):
    return math.sqrt(2 * math.log(12 / 0.01) / shots)


def graph_and_cells(board):
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    dist = nx.single_source_shortest_path_length(graph, 0)
    farthest = min(cell for cell, steps in dist.items() if steps == max(dist.values()))
    return graph, [0, farthest]


def source_and_reader(board, owner):
    block = sp.eye(12)
    for a, b in CELL_SEAMS:
        block[a, b] = block[b, a] = 1
    gates = [j for j in range(12) if board.degree[12 * owner + j] == 6]
    compact = block.col_join(sp.eye(12)[gates, :]).nullspace()[0]
    nums = [int(v * sp.ilcm(*(t.q for t in compact))) for v in compact]
    amps = np.zeros(board.n)
    amps[12 * owner : 12 * (owner + 1)] = nums
    amps /= np.linalg.norm(amps)
    support = np.flatnonzero(amps)
    choices = [
        (-int(nums[p % 12] ** 2), int(p), q, r)
        for p in support
        for q in board.neighbours[p]
        if not amps[q]
        for r in board.neighbours[q]
        if r != p and board.degree[r] == 5 and not amps[r]
    ]
    _, p, q, r = min(choices)
    return amps, {
        "owner": owner,
        "gateway_ports": gates,
        "integer_amplitudes_in_cell": nums,
        "sample_port_and_probe_rails": [p, q, r],
        "sample_vacancy_probability_if_in_cell": float(amps[p] ** 2),
    }


@cache
def signed_port_probabilities(a, b):
    column = law.echo()[:, law.SRC.index((a, b))]
    probabilities = []
    for d in probe.BINS:
        mask = np.repeat([False, *(v == d for _, v in law.LABELS)], law.N_SRC)
        probabilities.append(float(np.sum(abs(column[mask]) ** 2)))
    return np.asarray(probabilities)


def sharp_contrast(a, b):
    probabilities = signed_port_probabilities(a, b)
    maximum = int(np.argmax(probabilities))
    if abs(probabilities[maximum] - 1) > 1e-14:
        raise ValueError("A cut read has a non-deterministic outcome")
    return int(probe.BINS[maximum])


def cut_histories(board, graph, owner, first, second):
    reads = []
    for radius in (0, 1):
        inside = set(nx.single_source_shortest_path_length(graph, owner, cutoff=radius))
        edges = [
            (a, b) if a // 12 in inside else (b, a)
            for a, b in board.external
            if (a // 12 in inside) != (b // 12 in inside)
        ]
        endpoints = np.asarray(edges).ravel()
        if np.any(first[endpoints]) or np.any(second[endpoints]):
            raise ValueError("A charge amplitude reaches the measurement support")
        host_word = [sharp_contrast(1, 1) for _ in edges]
        challenge_vacancy = edges[0][0]
        challenge_word = [
            sharp_contrast(int(a != challenge_vacancy), int(b != challenge_vacancy))
            for a, b in edges
        ]
        reads.append(
            {
                "owner_radius": radius,
                "cut_oriented_from_inside": edges,
                "host_first_second_shared_signed_word": host_word,
                "vacancy_at_first_inner_endpoint_signed_word": challenge_word,
                "copies_each": SHOTS,
            }
        )
    return reads


def joint_histories(stations, curve, rng, shots=SHOTS, positive_only=False):
    curves = []
    for station in stations:
        weight = station["sample_vacancy_probability_if_in_cell"]
        curves.append((1 - weight) * curve[:, 1] + weight * curve[:, 0])
    histories = []
    for i, a in enumerate(probe.SIGNS):
        for j, b in enumerate(probe.SIGNS):
            if positive_only and (a, b) != (1, 1):
                continue
            host1, host2 = curve[i, 1], curve[j, 1]
            near1, near2 = curves[0][i], curves[1][j]
            shared = (np.outer(near1, host2) + np.outer(host1, near2)) / 2
            independent = np.outer((host1 + near1) / 2, (host2 + near2) / 2)
            for axis in (0, 1):
                np.testing.assert_allclose(
                    shared.sum(axis), independent.sum(axis), atol=1e-14
                )
            entry = {
                "probe_signs": [a, b],
                "joint_total_variation": float(abs(shared - independent).sum() / 2),
            }
            for label, probs in (
                ("shared", shared),
                ("independent_half_mean", independent),
            ):
                counts = probe.sample(shots, probs, rng)
                entry[label] = {
                    "counts": counts.tolist(),
                    "mean_score_product": float(SCORE @ counts @ SCORE / shots),
                }
            delta = (
                entry["independent_half_mean"]["mean_score_product"]
                - entry["shared"]["mean_score_product"]
            )
            entry["absolute_score_difference"] = abs(delta)
            histories.append(entry)
    return histories


def source_read(station, rng):
    """Signed echo of the sampled mode port against its first empty-of-defect rail."""
    w = station["sample_vacancy_probability_if_in_cell"]
    reads = {}
    for label, host, amount in (("vacancy", 1, 0), ("single", 0, 1), ("double", 0, 2)):
        probs = (1 - w) * signed_port_probabilities(
            host, host
        ) + w * signed_port_probabilities(amount, host)
        reads[label] = probe.sample(SHOTS, probs, rng).tolist()
    return reads


def conversion():
    mix = np.zeros((9, 9))
    pair = law.SRC.index((1, 1))
    for state in ((2, 0), (0, 2)):
        pos = law.SRC.index(state)
        mix[pos, pair] = mix[pair, pos] = K / np.sqrt(2)
    initial = law.SRC.index((2, 0))
    records = []
    rng = np.random.default_rng(SEED + 77)
    for value in (0, 1):
        column = expm(-1j / 16 * (free_seam(K) + value * mix))[:, initial]
        prob = [
            sum(abs(column[i]) ** 2 for i, (a, b) in enumerate(law.SRC) if a - b == d)
            for d in probe.BINS
        ]
        records.append(
            {
                "conversion_strength": float(value * K),
                "signed_word_counts": probe.sample(SHOTS, prob, rng).tolist(),
                "split_probability": float(prob[2]),
            }
        )
    return records


def run():
    rng = np.random.default_rng(SEED)
    precision_rng = np.random.default_rng(SEED + 1)
    cards = []
    curves = {k: probe.instrument(1 / 32, k) for k in (K, 0)}
    for level in range(3):
        board = federation(level)
        graph, owners = graph_and_cells(board)
        prepared = [source_and_reader(board, owner) for owner in owners]
        first, second = (entry[0] for entry in prepared)
        stations = [entry[1] for entry in prepared]
        for s in stations:
            rails = s["sample_port_and_probe_rails"][1:]
            assert not np.any(first[rails]) and not np.any(second[rails])
        shared = (first + second) / np.sqrt(2)
        eigen = {}
        for name, h in (
            ("vacancy_in_host", defect_generator(board, -1, K, 0)),
            (
                "single_on_empty",
                board.laplacian() / 2 + 6 * K * eye(board.n, format="csr"),
            ),
            ("double_on_empty", defect_generator(board, 2, K, 0, background=0)),
        ):
            e = float(shared @ (h @ shared))
            residual = max(
                float(np.max(abs(h @ v - e * v))) for v in (first, second, shared)
            )
            if residual > 2e-12:
                raise ValueError("Source leaks during free release")
            eigen[name] = {
                "energy_above_host": e,
                "largest_component_residual": residual,
            }

        for owner, s in zip(owners, stations, strict=True):
            s["cuts"] = cut_histories(board, graph, owner, first, second)
        cards.append(
            {
                "cells": board.cells,
                "owners_apart": nx.shortest_path_length(graph, *owners),
                "stations": stations,
                "free_eigenstates": eigen,
                "source_signed_read": source_read(
                    stations[0], rng=np.random.default_rng(SEED + 10 + level)
                ),
                "precision_positive_probes": joint_histories(
                    stations,
                    curves[K],
                    precision_rng,
                    shots=PRECISION_SHOTS,
                    positive_only=True,
                )[0],
                "joint_probe_settings": [
                    {
                        "terminal_contact_strength": strength,
                        "histories": joint_histories(stations, curve, rng),
                    }
                    for strength, curve in curves.items()
                ],
            }
        )
    return {
        "shots_per_setting": SHOTS,
        "seed": SEED,
        "one_mean_simultaneous_99_percent_error": simultaneous_radius(SHOTS),
        "precision_shots": PRECISION_SHOTS,
        "precision_seed": SEED + 1,
        "precision_one_mean_simultaneous_99_percent_error": simultaneous_radius(
            PRECISION_SHOTS
        ),
        "positive_probe_calibration": {
            "signed_contrast_bins": probe.BINS.tolist(),
            "ternary_input_probabilities": curves[K][1].tolist(),
        },
        "isolated_symmetric_total_two_challenge": conversion(),
        "boards": cards,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = json.loads(json.dumps(run()))
    if args.output:
        with args.output.open("x") as f:
            json.dump(receipt, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(receipt, json.loads((HERE / "results.json").read_text()), "split")
        print(
            "Reproduced the two contact histories, enclosing cuts and conversion control."
        )
