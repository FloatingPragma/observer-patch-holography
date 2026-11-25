"""Can the entire written-clock response be explained by arbitrary random elapsed times?"""

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_03_clock.written_coin import measure as clock
from model.candidate_laws.interval_echo.openwave.row_27_metric.free_source_clock import measure as free
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label

HERE = Path(__file__).resolve().parent
A, B = free.ENDPOINTS
WAITS = np.array([0, .25, .5, 1, 2, 4, 8, 12])
SHOTS, TIMER_SHOTS = clock.SHOTS, clock.TIMER_SHOTS
SOURCE_SEED = 27112510
# Two laws, two orientations, eleven source settings, five material bins plus timer.
FAMILY = 2 * 2 * 11 * len(WAITS) * 6
TAIL = .01 / (2 * FAMILY)
EPS = np.sqrt(-np.log(TAIL) / (2 * SHOTS))
READS = [(2, 0), (0, 2), (2, 1), (1, 2)]


def local_blocks(control):
    law = candidate.native_law
    e = law.echo()
    exchange, contrast = law.generator_parts()
    if control:
        magnitude_two = [1 + i for i, (_, d) in enumerate(law.LABELS) if abs(d) == 2]
        removed = np.concatenate([i * law.N_SRC + np.arange(law.N_SRC) for i in magnitude_two])
        exchange[removed, :] = 0
        exchange[:, removed] = 0
    h = e @ (exchange + candidate.KAPPA * contrast) @ e.T
    if not control:
        ex, co = law.local_term_parts()
        np.testing.assert_allclose(h, ex + candidate.KAPPA * co, rtol=0, atol=1e-13)
    loads = np.eye(law.N_SRC)
    labels = [label(2, p, 1) for p in (1, -1)]
    gaps = (4 * candidate.KAPPA, 4 * candidate.KAPPA + (0 if control else 1))
    for sign, orientation in zip((1, -1), READS[:2], strict=True):
        memory = (labels[0] + sign * labels[1]) / np.sqrt(2)
        initial = np.eye(law.DIM)[:, law.SRC.index(orientation)]
        for local in ((1, 1), *READS[2:]):
            v = np.kron(memory, loads[law.SRC.index(local)])
            if local == (1, 1):
                np.testing.assert_allclose(e @ initial, v, atol=1e-13)
                np.testing.assert_allclose(e @ v, initial, atol=1e-13)
            else:
                np.testing.assert_allclose(e @ v, v, atol=1e-13)
            for arm, gap in zip(labels, gaps, strict=True):
                lifted = np.kron(arm, loads[law.SRC.index(local)])
                np.testing.assert_allclose(h @ lifted, (0 if local == (1, 1) else gap) * lifted,
                                           rtol=0, atol=1e-13)
    # Every other seam sees either filled endpoints or a single excess. The
    # control preserves all columns needed for its free propagation, including leakage.
    for local in ((1, 1), *READS[2:]):
        i = law.SRC.index(local)
        np.testing.assert_allclose(h[:, i], np.r_[law.seam_term(candidate.KAPPA)[:, i],
                                                  np.zeros(law.DIM - law.N_SRC)], atol=1e-13)
    return gaps


def branch(h, starts):
    energies, v = eigh(h, driver="evd")
    phase = np.exp(-1j * WAITS[:, None] * energies)
    return {j: (phase * v[j, :]) @ v.T for j in starts}


def bins(even, odd):
    same, reverse = (abs(even + odd) ** 2 / 4, abs(even - odd) ** 2 / 4)
    away = np.ones(even.shape[1], dtype=bool)
    away[[A, B]] = False
    probabilities = np.column_stack((same[:, away].sum(1), reverse[:, away].sum(1),
                                    same[:, A] + reverse[:, A], same[:, B] + reverse[:, B],
                                    np.zeros(len(WAITS))))
    np.testing.assert_allclose(probabilities.sum(axis=1), 1, rtol=0, atol=3e-13)
    return probabilities


def interval(count):
    lo = float(beta.ppf(TAIL, count, SHOTS - count + 1)) if count else 0.
    hi = float(beta.isf(TAIL, count + 1, SHOTS - count)) if count < SHOTS else 1.
    return lo, hi


def return_certificate(adjacency, site, distance):
    """Exact leading opposite-sign probability with the excess ending at launch.

    At an endpoint use a specific adjacent away port instead, since the local
    endpoint read cannot resolve sign while the excess is there.
    """
    path_counts = np.eye(1, len(adjacency), site, dtype=object).ravel()
    for _ in range(distance):
        path_counts = adjacency @ path_counts
    if distance:
        overlap = sum(int(path_counts[v]) ** 2 for v in (A, B))
        order, hop_weight, terminal = 2 * distance + 1, Fraction(overlap, 4**distance), site
    else:
        terminal = next(i for i in np.flatnonzero(adjacency[site]) if i not in (A, B))
        order, hop_weight = 2, Fraction(1, 2)
    coeff = (hop_weight / (2 * factorial(order))) ** 2
    if not coeff:
        raise ValueError("missing shortest-path amplitude")
    return {"terminal_excess_port": int(terminal), "probability_power": 2 * order,
            "leading_probability_coefficient": str(coeff)}


def acquire():
    board = federation(0)
    adjacency = board.adjacency().toarray().astype(int)
    dist = shortest_path(adjacency, directed=False, indices=[A, B], unweighted=True).min(0)
    sites = {d: int(np.flatnonzero(dist == d)[0]) for d in (0, 1, 2, 3, int(max(dist)))}
    adjacency[A, B] = adjacency[B, A] = 0
    h = free.board_generator(board, candidate.KAPPA)
    mask = np.diag([int(i in (A, B)) for i in range(board.n)])
    rng = np.random.default_rng(SOURCE_SEED)
    raw, summary = {}, {}
    for name, control in (("native", False), ("no_sign_exchange_at_two", True)):
        gaps = local_blocks(control)
        first = branch(h + gaps[0] * mask, sites.values())
        second = branch(h + gaps[1] * mask, sites.values()) if gaps[1] != gaps[0] else first
        probs = {"blank": np.tile([1., 0, 0, 0, 0], (len(WAITS), 1))}
        for d, j in sites.items():
            probs[f"free{d}"] = bins(first[j], second[j])
            # Shuttered ordinary links hold the excess fixed. Its global
            # phase is immaterial; an occupied endpoint gives an unresolved read.
            parked = np.zeros(5)
            parked[{A: 2, B: 3}.get(j, 0)] = 1.
            probs[f"park{d}"] = np.tile(parked, (len(WAITS), 1))
        for sign in ("plus", "minus"):
            for arm, p in probs.items():
                key = f"{name}_{sign}_{arm}"
                physical = p if sign == "plus" else p[:, [1, 0, 2, 3, 4]]
                bank = np.array([clock.enroll(rng, row, SHOTS) for row in physical])
                raw[key] = bank
                raw[key + "_timer"] = rng.binomial(TIMER_SHOTS, np.sin(WAITS / 16) ** 2)
                reverse = bank[:, 1 if sign == "plus" else 0]
                blank = raw[f"{name}_{sign}_blank"][:, 1 if sign == "plus" else 0]
                resolved = [interval(r)[0] > interval(b)[1] for r, b in zip(reverse, blank, strict=True)]
                summary[key] = {"resolved_waits_exact_binomial": WAITS[resolved].tolist(),
                                "resolved_waits_hoeffding": WAITS[(reverse - blank) / SHOTS > 2 * EPS].tolist()}
    paths = {str(d): return_certificate(adjacency, j, d) for d, j in sites.items()}
    return ({"endpoints": [A, B], "launch_by_native_distance": sites, "waits": WAITS.tolist(),
             "shots": SHOTS, "timer_shots": TIMER_SHOTS, "count_columns": [*map(list, READS), "other"],
             "simultaneous_99_probability_radius": float(EPS), "comparison_marginals": FAMILY,
             "held_free1_waits": WAITS[3:].tolist(), "short_time_certificates": paths,
             "decisions": summary}, raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    summary, raw = acquire()
    report = json.dumps(summary, indent=2) + "\n"
    if args.write:
        np.savez_compressed(HERE / "counts.npz", **raw)
        (HERE / "results.json").write_text(report)
    else:
        expected = json.loads((HERE / "results.json").read_text())
        if json.loads(report) != expected:
            raise ValueError("summary differs")
        with np.load(HERE / "counts.npz") as stored:
            if set(raw) != set(stored.files):
                raise ValueError("record inventory differs")
            for key, value in raw.items():
                np.testing.assert_array_equal(value, stored[key])
    for sign, col in (("plus", 1), ("minus", 0)):
        key = f"native_{sign}_free1"
        print(key, raw[key][3:, col].tolist(), summary["decisions"][key])


if __name__ == "__main__":
    main()
