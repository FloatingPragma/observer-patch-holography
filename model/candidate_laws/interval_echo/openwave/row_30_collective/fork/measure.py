"""Simulate free strand-wall scattering on three physical branch directions."""

import argparse
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.optimize import brentq
from scipy.stats import binom

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30_collective.experiment import (
    native_checks,
)
from model.executable_model import federation

HERE = Path(__file__).parent
EXPOSURES = (0, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0)
SHOTS = 2**28
TAIL = 0.01 / (
    2 * 50000
)  # covers every reported event/comparison, including phase settings.


def fixture(board, length):
    center = (6, 7)
    across = dict(list(board.external) + [(b, a) for a, b in board.external])
    gates = sorted(g for g in across if g // 12 == 0)
    dimers, ports, links = [center], [], []
    used = set(center)
    for gate in gates:
        foot = min(p for p in center if p in board.neighbours[gate])
        previous, tip = 0, foot
        for j in range(length):
            if j == 0:
                pair = (gate, across[gate])
            else:
                a = min(
                    x
                    for x in board.neighbours[tip]
                    if x not in used and x not in across and x // 12 == tip // 12
                )
                b = min(
                    x
                    for x in board.neighbours[a]
                    if x not in used and x not in across and x // 12 == tip // 12
                )
                pair = (a, b)
            if (
                any(p in used for p in pair)
                or board.degree[pair[0]] != board.degree[pair[1]]
            ):
                raise ValueError("dimers must be disjoint and degree matched")
            if (
                pair[1] not in board.neighbours[pair[0]]
                or pair[0] not in board.neighbours[tip]
            ):
                raise ValueError("non-native edge")
            i = len(dimers)
            dimers.append(pair)
            ports.append((tip, pair[0]))
            links.append((previous, dimers[previous].index(tip), i, 0))
            used.update(pair)
            previous, tip = i, pair[1]
    if len({pair[1] // 12 for pair in dimers[1::length]}) != 3:
        raise ValueError("the three gateways must reach distinct receiver owners")
    return dimers, ports, links


def material(links, kappa):
    n = len(links) + 1
    words = np.arange(2**n)
    orientations = words[:, None] >> np.arange(n) & 1
    walls = np.zeros_like(words)
    cost = np.zeros_like(words)
    for j, (a, pa, b, pb) in enumerate(links):
        left, right = orientations[:, a] == pa, orientations[:, b] == pb
        walls += (left != right).astype(int) << j
        cost += 1 - (left & right)
    # Matched degrees make the private and internal reference energies scalars.
    # The matrix omits the common diagonal of the dimer exchange and collar terms.
    rows = np.tile(words, n)
    cols = np.concatenate([words ^ (1 << i) for i in range(n)])
    h = sparse.csr_matrix(
        (np.full(rows.size, -0.5), (rows, cols)), shape=(len(words), len(words))
    )
    h += sparse.diags(4 * kappa * (cost - cost.min()), dtype=float)
    return h, walls


def events(length):
    return {
        "input": 1,
        "along_A": 2,
        "turn_B": 1 << length,
        "turn_C": 1 << (2 * length),
        "split_at_center": (1 << length) | (1 << (2 * length)),
        "split_in_B_and_C": (2 << length) | (2 << (2 * length)),
    }


def evolve(h, launch):
    initial = np.zeros(h.shape[0], dtype=complex)
    initial[launch] = 1.0
    result = [initial]
    for before, after in zip(EXPOSURES, EXPOSURES[1:]):
        result.append(expm_multiply(-1j * (after - before) * h, result[-1]))
    return np.array(result)


def acquire(psi, rng):
    p = np.round(np.abs(psi) ** 2, 13)
    p /= p.sum()
    return rng.multinomial(SHOTS, p)


@lru_cache(maxsize=8192)
def band(k):
    # Invert both tails directly: some extreme beta quantiles at this shot
    # budget can fall on the wrong side of the observed proportion.
    if k > SHOTS // 2:
        low, high = band(SHOTS - k)
        return [
            max(0.0, float(np.nextafter(1 - high, 0.0))),
            min(1.0, float(np.nextafter(1 - low, 1.0))),
        ]
    center = k / SHOTS
    lo = (
        brentq(lambda p: binom.sf(k - 1, SHOTS, p) - TAIL, 0, center, xtol=1e-24)
        if k
        else 0.0
    )
    hi = (
        brentq(lambda p: binom.cdf(k, SHOTS, p) - TAIL, center, 1, xtol=1e-24)
        if k != SHOTS
        else 1.0
    )
    # Round away from the measured proportion, including solver tolerance.
    margin = 1e-24 + 2e-15 * hi
    return [max(0.0, lo - margin), min(1.0, hi + margin)]


def paired_interval(a, b):
    aa, bb = band(a), band(b)
    return [aa[0] - bb[1], aa[1] - bb[0]]


def interval_distance(interval):
    return max(0.0, interval[0], -interval[1])


def run():
    native = native_checks()
    rng = np.random.default_rng(790301)
    counts, summary = (
        {},
        {
            "native_column_residual": native["max_amplitude_residual"],
            "shots": SHOTS,
            "tail_per_marginal": TAIL,
            "fixtures": {},
            "runs": {},
        },
    )
    for level in range(3):
        board = federation(level)
        for length in (2, 3):
            dimers, ports, links = fixture(board, length)
            tag = f"b{board.cells}_L{length}"
            summary["fixtures"][tag] = {"dimers": dimers, "connectors": ports}
            for k in (candidate.KAPPA, 0):
                h, walls = material(links, k)
                empty = int(next(x for x in np.flatnonzero(walls == 0) if not x & 1))
                incident = empty ^ sum(1 << d for d in range(1, length + 1))
                assert walls[incident] == 1
                masks = events(length)
                for reverse in (0, 1):
                    key = f"{tag}_k{int(k)}_r{reverse}"
                    both = (len(walls) - 1) * reverse
                    panel = {}
                    for name, launch in [
                        ("incident", incident ^ both),
                        ("blank", empty ^ both),
                    ]:
                        flight = evolve(h, launch)
                        table = np.empty(
                            (len(EXPOSURES), 3, len(walls)), dtype=np.int64
                        )
                        for t, psi in enumerate(flight):
                            np.testing.assert_allclose(
                                np.vdot(psi, psi), 1.0, atol=2e-10
                            )
                            other = psi[np.arange(len(walls)) ^ 1]
                            for arm, state in enumerate(
                                [
                                    psi,
                                    (psi + 1j * other) / np.sqrt(2),
                                    (psi - 1j * other) / np.sqrt(2),
                                ]
                            ):
                                table[t, arm] = acquire(state, rng)
                        counts[f"{key}_{name}"] = table
                        panel[name] = {
                            m: table[:, 0, :][:, walls == mask].sum(1)
                            for m, mask in masks.items()
                        }
                        remote = masks["split_in_B_and_C"]
                        panel[name]["remote_both"] = table[:, 0, :][
                            :, (walls & remote) == remote
                        ].sum(1)
                        panel[name]["phase"] = table[:, 1:, launch]
                    timer = rng.binomial(SHOTS, np.sin(np.array(EXPOSURES) / 8) ** 2)
                    counts[key + "_timer"] = timer
                    result = {}
                    for ev in (*masks, "remote_both"):
                        a, b = panel["incident"][ev], panel["blank"][ev]
                        bounds = [paired_interval(int(n), int(m)) for n, m in zip(a, b)]
                        i = max(
                            range(len(bounds)),
                            key=lambda j: interval_distance(bounds[j]),
                        )
                        result[ev] = dict(
                            index=i,
                            incident=int(a[i]),
                            blank=int(b[i]),
                            difference_band=bounds[i],
                            all_incident=a.tolist(),
                            all_blank=b.tolist(),
                        )
                    plus, minus = panel["incident"]["phase"].T
                    bounds = [
                        paired_interval(int(n), int(m)) for n, m in zip(plus, minus)
                    ]
                    i = max(
                        range(len(bounds)), key=lambda j: interval_distance(bounds[j])
                    )
                    result["center_coherence"] = dict(
                        index=i,
                        plus=int(plus[i]),
                        minus=int(minus[i]),
                        difference_band=bounds[i],
                    )
                    result["times_from_hand"] = (
                        8 * np.arcsin(np.sqrt(timer / SHOTS))
                    ).tolist()
                    summary["runs"][key] = result
    return json.loads(json.dumps(summary)), counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result, counts = run()
    if args.output:
        args.output.mkdir(exist_ok=False)
        (args.output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
        np.savez_compressed(args.output / "counts.npz", **counts)
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(result, json.loads((HERE / "results.json").read_text()), "fork")
        with np.load(HERE / "counts.npz") as stored:
            assert set(counts) == set(stored.files)
            for key, bank in counts.items():
                np.testing.assert_array_equal(bank, stored[key])
        print("Replayed every fork orientation, phase and clock bank.")


if __name__ == "__main__":
    main()
