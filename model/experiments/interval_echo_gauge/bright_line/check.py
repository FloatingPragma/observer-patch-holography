"""Check all records using the raw permutation Hamiltonian and a spectral projection."""

from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

from model.executable_model import CELL_SEAMS


def direct(loads, receipt):
    edges = [tuple(sorted(e)) for e in sorted(CELL_SEAMS)]
    adjacency = np.zeros((12, 12), dtype=int)
    for a, b in edges:
        adjacency[a, b] = adjacency[b, a] = 1
    anti = shortest_path(adjacency, directed=False).argmax(axis=1)
    words = list(combinations(range(12), loads))
    h = np.zeros((len(words), len(words)))
    sources = np.zeros((len(words), 12))
    for k, word in enumerate(words):
        sites = set(word)
        for a, b in edges:
            if {a, b} <= sites:
                h[k, k] -= 16
            if (a in sites) != (b in sites):
                other = tuple(sorted(sites ^ {a, b}))
                h[k, k] += 0.5
                h[words.index(other), k] -= 0.5
        for p in range(12):

            def launch(q):
                return (q in sites) and (loads == 1 or word in edges)

            sources[k, p] = int(launch(p)) - int(launch(anti[p]))
    sources /= np.linalg.norm(sources, axis=0)
    eig, columns = np.linalg.eigh(h)
    active = np.where(abs(columns.T @ sources[:, 1]) > 1e-7)[0]
    first = columns[:, abs(eig - eig[active[0]]) < 1e-9]
    assert first.shape[1] == 3
    occupation = np.array([len(set(word) & {1, anti[1]}) for word in words])
    d = first.T @ (occupation[:, None] * first)
    vals = np.linalg.eigvalsh(d)
    np.testing.assert_allclose(vals[0], vals[1], atol=1e-12)
    beta = vals[2] if loads == 2 else sum(receipt["pair_bright_enclosure"]) / 2
    times = dict(
        sham=(0, 0),
        phase_anchor=(0, np.pi / loads),
        native=(np.pi, 0),
        determinant_lock=(np.pi, np.pi * (beta - 0.5)),
        fixed_line=(np.pi, -np.pi * beta / 2),
        plane_flip=(4 * np.pi * receipt["flip_cycles"], 0),
        zero_occupation=(0, 0),
    )
    if loads == 2:
        ordered = np.array([words.index(edge) for edge in edges])
        quiet = np.ones(len(words), dtype=bool)
        quiet[ordered] = False
    enroll = receipt["shots"]
    output = receipt["single" if loads == 1 else "pair"]
    for arm, (time, cal) in times.items():
        u = expm(1j * (time * d + loads * cal * np.eye(3)))
        if arm == "plane_flip":
            bright_vector = first.T @ sources[:, 1]
            bright_vector /= np.linalg.norm(bright_vector)
            target = (
                np.eye(3)
                if loads == 1
                else 2 * np.outer(bright_vector, bright_vector) - np.eye(3)
            )
            err = np.linalg.norm(u - target, 2)
            assert err < (1e-9 if loads == 1 else 0.033)
        final = sources + first @ ((u - np.eye(3)) @ (first.T @ sources))
        for k in (0, 1):
            prob = abs(final[:, k]) ** 2
            if loads == 2:
                prob = np.r_[prob[ordered], sum(prob[quiet])]
            data = output[arm]["launches"][str(k)]
            # Finite symmetric sum-pulse errors plus printed calibration enclosure.
            np.testing.assert_allclose(data["probabilities"], prob, rtol=0, atol=7e-6)
            counts = np.array(data["counts"])
            assert counts.sum() == enroll and counts.min() >= 0
            assert (
                np.max(abs(counts / enroll - prob))
                < np.sqrt(np.log(2 * (7 * 2 * (12 + 31)) / 0.01) / (2 * enroll)) + 7e-6
            )


if __name__ == "__main__":
    receipt = json.loads(Path(__file__).with_name("records.json").read_text())
    for loads in (1, 2):
        direct(loads, receipt)
    print(
        "All 28 full observer distributions agree with the independent raw-cell calculation."
    )
