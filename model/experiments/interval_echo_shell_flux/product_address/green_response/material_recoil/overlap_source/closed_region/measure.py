"""Closed native energy and occupancy after two separated records are prepared."""

import argparse
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import CELL_SEAMS, federation

HERE = Path(__file__).parent
TIMES = (0, 1, 4, 16, 64, 256)


def layout(owners):
    if owners == 1:
        ports, edges, anchor = list(range(12)), CELL_SEAMS, 0
    else:
        board = federation(1)
        a, b = min(board.external)
        ports = [12 * owner + j for owner in (a // 12, b // 12) for j in range(12)]
        pos = {v: i for i, v in enumerate(ports)}
        edges = [(pos[a], pos[b]) for a, b in board.seams if a in pos and b in pos]
        anchor = pos[a]
    graph = np.zeros((len(ports), len(ports)), int)
    for a, b in edges:
        graph[a, b] = graph[b, a] = 1
    distance = shortest_path(graph, directed=False, unweighted=True)[anchor].astype(int)
    adjacent = [i for i in range(12) if graph[anchor, i]]
    starts = [pair for pair in combinations(adjacent, 2) if not graph[pair]]
    assert len(starts) == 5
    return ports, edges, distance, starts


@dataclass
class Sector:
    words: list
    matrix: np.ndarray  # twice the complete native Hamiltonian, including its diagonal
    indices: list
    charges: np.ndarray
    cost: int


def sector(n, edges, count, k):
    words = list(combinations(range(n), count))
    masks = [sum(1 << v for v in word) for word in words]
    reverse = {mask: j for j, mask in enumerate(masks)}
    mat = np.zeros((len(words), len(words)), dtype=np.int64)
    indices = []
    for a, b in edges:
        left = np.array(
            [i for i, m in enumerate(masks) if m & (1 << a) and not m & (1 << b)], int
        )
        right = np.array([reverse[masks[i] ^ (1 << a) ^ (1 << b)] for i in left], int)
        mat[left, left] += 1 + 2 * k
        mat[right, right] += 1 + 2 * k
        mat[left, right] -= 1
        mat[right, left] -= 1
        indices.append((left, right))
    charges = np.array([[int(v in word) for word in words] for v in range(n)])
    return Sector(words, mat, indices, charges, 1 + 2 * k)


def observables(sec, rho):
    diagonal = np.diag(rho).real
    energy = np.array(
        [
            (sec.cost * (diagonal[i] + diagonal[j]).sum() / 2 - rho[i, j].real.sum())
            for i, j in sec.indices
        ]
    )
    return energy, sec.charges @ diagonal


def native_check(k):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    ix = [law.SRC.index(w) for w in ((0, 0), (0, 1), (1, 0), (1, 1))]
    embed = np.eye(law.DIM)[:, ix]
    expected = np.zeros((4, 4))
    expected[1:3, 1:3] = [[1 + 2 * k, -1], [-1, 1 + 2 * k]]
    np.testing.assert_allclose(
        2 * (exchange + k * contrast) @ embed, embed @ expected, atol=4e-14
    )


def readings(energy, charge, edges, distance):
    shells = []
    for r in range(int(distance.max()) + 1):
        inside = distance <= r
        share = np.array([(int(inside[a]) + int(inside[b])) / 2 for a, b in edges])
        shells.append(
            dict(
                radius=r, energy=float(share @ energy), load=float(charge[inside].sum())
            )
        )
    return shells


def evolved(sec, cohorts):
    eigen, vect = eigh(sec.matrix / 2)
    # Exact eigenvalue multiplicities on the cell are checked independently.
    # The two-cell asymptotic forecast is numerical; expose its grouping sensitivity.
    energies = {}
    diagnostics = {}
    for label, starts in cohorts.items():
        initial = np.zeros(len(sec.words))
        for i, weight in starts:
            initial[i] += weight
        base = (vect.T * initial) @ vect
        rows = {}
        for t in TIMES:
            phase = np.exp(-1j * t * eigen)
            rho = vect @ (phase[:, None] * base * phase.conj()[None, :]) @ vect.T
            rows[str(t)] = observables(sec, rho)
        gaps = eigen[:, None] - eigen[None, :]
        for duration in (16, 256, 65536):
            angle = gaps * duration / 2
            window = np.exp(-1j * angle) * np.sinc(angle / np.pi)
            rho = vect @ (window * base) @ vect.T
            rows[f"mean_through_{duration}"] = observables(sec, rho)
        for tol in (1e-11, 1e-9):
            groups = np.cumsum(np.r_[False, np.diff(eigen) > tol])
            block = base * (groups[:, None] == groups[None, :])
            answer = observables(sec, vect @ block @ vect.T)
            if tol == 1e-11:
                rows["mean"] = answer
            else:
                diagnostics[label] = max(
                    float(abs(a - b).max())
                    for a, b in zip(rows["mean"], answer, strict=True)
                )
        energies[label] = rows
    return energies, diagnostics


def run_case(owners, k):
    native_check(k)
    ports, edges, distance, pairs = layout(owners)
    singles = sector(len(ports), edges, 1, k)
    doubles = sector(len(ports), edges, 2, k)
    names = [str((ports[a], ports[b])) for a, b in pairs]
    s = {
        key: [(singles.words.index((v,)), 1) for v in pair]
        for key, pair in zip(names, pairs, strict=True)
    }
    t = {
        key: [(doubles.words.index(pair), 1)]
        for key, pair in zip(names, pairs, strict=True)
    }
    a, sa = evolved(singles, s)
    b, sb = evolved(doubles, t)
    rows = {}
    for key in names:
        rows[key] = {}
        for time, one in a[key].items():
            en, occ = (
                single - double
                for single, double in zip(one, b[key][time], strict=True)
            )
            assert abs(en.sum()) < 3e-11 and abs(occ.sum()) < 3e-12
            if time == "0":
                assert max(abs(en)) < 3e-11
            rows[key][time] = readings(en, occ, edges, distance)
    average = {}
    ordinary = {}
    for time in a[names[0]]:
        # One ordinary word, with its starting endpoint drawn uniformly from
        # the two input sites and all enrolled placements. The vacuum gives zero.
        ordinary_energy, ordinary_charge = (
            np.mean([a[key][time][0] / 2 for key in names], axis=0),
            np.mean([a[key][time][1] / 2 for key in names], axis=0),
        )
        ordinary[time] = readings(ordinary_energy, ordinary_charge, edges, distance)
        profiles = [rows[key][time] for key in names]
        average[time] = [
            dict(
                radius=i,
                energy=float(np.mean([p[i]["energy"] for p in profiles])),
                load=float(np.mean([p[i]["load"] for p in profiles])),
            )
            for i in range(int(distance.max()) + 1)
        ]
    return dict(
        cells=owners,
        stiffness=k,
        native_ports=ports,
        pairs=names,
        energy_unit="native energy per unit of original covariance",
        maximum_degeneracy_threshold_effect=max((*sa.values(), *sb.values())),
        average=average,
        ordinary_single_load=ordinary,
        individual=rows,
    )


def run():
    inputs = json.loads(HERE.parent.joinpath("formed_heads.json").read_text())
    return dict(
        covariance={
            str(x["source_sign"]): x["covariance_of_original_load_bits"]
            for x in inputs["source"]
            if x["orientation"] == "matching"
        },
        panels=[run_case(n, k) for n in (1, 2) for k in (8, 0)],
    )


def check(old, new):
    if isinstance(old, dict):
        assert old.keys() == new.keys()
        for key in old:
            check(old[key], new[key])
    elif isinstance(old, list):
        for a, b in zip(old, new, strict=True):
            check(a, b)
    elif isinstance(old, float):
        np.testing.assert_allclose(old, new, rtol=1e-8, atol=2e-10)
    else:
        assert old == new


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    path = HERE / "closed.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), data)
