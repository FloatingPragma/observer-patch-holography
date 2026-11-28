"""Measure all station times and both ordinary energies before inferring a regional source."""

import argparse
from hashlib import sha256
from itertools import combinations
import json

import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import paired_balls as balls

survey = balls.survey
HERE = balls.HERE


def partition(records, ref, independent, delta):
    _, labels = np.unique(np.ascontiguousarray(records), axis=0, return_inverse=True)

    def group(p):
        return np.bincount(labels, weights=p.ravel()).astype(int)

    p, f, d = map(group, (ref, independent, delta))

    def entropy(p):
        return float(-np.sum(p / 2592 * np.log(p / 2592)))

    return {
        "classes": len(p),
        "sum_absolute_count_change": int(abs(d).sum()),
        "maximum_discrimination_fraction": str(balls.F(int(abs(d).sum()), 72)),
        "entropy_slope": float(-(d / 2592) @ np.log(p / 2592)),
        "finite_entropy_change": entropy(f) - entropy(p),
    }, labels


def fixed_contact(calibration, price):
    """A preset seam pulse and terminal tag position, independent of the word."""
    from scipy.linalg import expm

    edges = sorted(survey.routes.cell().seams)
    bits = np.array(calibration["1"][1:])
    hits = [j for j in range(len(edges)) if len(set(zip(price, bits[:, j]))) == 9]
    head = hits[0]
    edge = edges[head]
    volume = survey.routes.volume
    law = volume.candidate.native_law
    exchange, cost = law.local_term_parts()
    h = exchange + volume.candidate.KAPPA * cost
    unitary = expm(-1j * np.pi * h)
    labels = [volume.label(0, 1), volume.label(2, 1), volume.label(2, -1)]
    observed, energies = [], []
    blank = np.eye(11)[0]
    vacuum = np.eye(law.N_SRC)[law.SRC.index((0, 0))]
    for k, matching in enumerate(survey.routes.codebook()[0][1:]):
        letters = divmod(k, 3)
        marks = {
            seam: labels[letters[n] if n < 2 else 0] for n, seam in enumerate(matching)
        }
        # A single fixed station-energy instrument sums every internal seam
        # on the parked preparation, including all blank contacts. Each
        # input is a simultaneous eigenpreparation before the tag launches.
        local_energy = 0
        for seam in edges:
            sharp = np.kron(marks.get(seam, blank), vacuum)
            applied = h @ sharp
            value = round(float(np.vdot(sharp, applied).real))
            np.testing.assert_allclose(applied, value * sharp, rtol=0, atol=3e-13)
            local_energy += value
        assert local_energy == price[k]
        energies.append(local_energy)
        memory = marks.get(edge, blank)
        row = []
        for q in (1, 2):
            init = np.kron(memory, np.eye(law.N_SRC)[law.SRC.index((q, 0))])
            out = (unitary @ init).reshape(11, law.N_SRC)
            arrive = float(np.linalg.norm(out[:, law.SRC.index((0, q))]) ** 2)
            np.testing.assert_allclose(arrive, bits[k, head], rtol=0, atol=3e-13)
            # The timing pulse preserves the sharp memory up to a phase.
            np.testing.assert_allclose(
                out, np.outer(memory, memory.conj() @ out), rtol=0, atol=3e-13
            )
            row.append(arrive)
        observed.append(row)
    return {
        "ports": [list(edge)],
        "minimal_decks": len(hits),
        "sharp_records_in_calibration_order": np.column_stack(
            [energies, np.rint(np.asarray(observed)[:, 0]).astype(int)]
        ).tolist(),
        "one_and_two_load_arrival_probabilities": observed,
    }


def local_reader(calibration, price, stocks, rng):
    """Calibrate the detector on sharp preparations before comparing source arms."""
    ports = list(combinations(range(12), 2))
    edges = sorted(survey.routes.cell().seams)
    times = []
    for outcomes in calibration["1"][1:]:
        graph = np.zeros((12, 12), dtype=int)
        for (i, j), passed in zip(edges, outcomes, strict=True):
            graph[i, j] = graph[j, i] = passed
        times.append(shortest_path(graph, unweighted=True)[np.triu_indices(12, 1)])
    times = np.array(times, dtype=int)
    decks = {}
    for kind, pool in (("all", range(len(ports))), ("port_zero", range(11))):
        for size in range(len(pool) + 1):
            hits = []
            for subset in combinations(pool, size):
                symbols = np.column_stack([price, times[:, subset]])
                if len(np.unique(symbols, axis=0)) == 9:
                    hits.append(subset)
            if hits:
                use = hits[0]
                symbols = np.column_stack([price, times[:, use]])
                decks[kind] = {
                    "ports": [list(ports[i]) for i in use],
                    "minimal_decks": len(hits),
                    "sharp_records_in_calibration_order": symbols.tolist(),
                }
                break
    decks["fixed_contact"] = fixed_contact(calibration, price)
    # Model the observed histogram directly in lexicographic instrument order.
    symbols = np.array(decks["fixed_contact"]["sharp_records_in_calibration_order"])
    pairs = np.array([np.r_[u, v] for u in symbols for v in symbols])
    alphabet, labels = np.unique(pairs, axis=0, return_inverse=True)
    assert len(alphabet) == 81
    n = 2**50
    iv = balls.IV
    iv.dps = 35
    trial = {}
    for tag, triple in stocks:
        counts = []
        for law in triple:
            numerators = np.bincount(labels, weights=law.ravel()).astype(int)
            counts.append(rng.multinomial(n, numerators / 2592).tolist())
        # Union bound over all 486 bin frequencies (both encodings).
        radius = iv.sqrt(iv.log(2 * 486 / iv.mpf("0.01")) / (2 * n))
        error = iv.mpf([-1, 1]) * radius
        rows = [[iv.mpf(c) / n + error for c in row] for row in counts]
        ref, fresh, shared = rows
        slope = -sum((v - u) * iv.log(81 * b) for b, u, v in zip(*rows, strict=True))

        def relative_entropy_term(p):
            return p * iv.log(81 * p)

        finite = sum(
            relative_entropy_term(b) - relative_entropy_term(u)
            for b, u in zip(ref, fresh, strict=True)
        )
        trial[tag] = {
            "reference_fresh_shared_counts": counts,
            "modular_slope": survey.routes.bounds(slope),
            "finite_entropy_increase": survey.routes.bounds(finite),
        }
    return {
        "calibrations": decks,
        "primary": "fixed_contact",
        "alphabet": alphabet.tolist(),
        "shots_per_arm": n,
        "simultaneous_coverage": 0.99,
        "acquisitions": trial,
    }


def run():
    board = survey.bulk.federation(0)
    weights, turned = survey.standards(board)
    src = survey.routes.volume.radar_source_and_hand(board)[-1]
    energy = np.array(src["source_word_energies"])
    price = (energy[:, None] + energy).ravel()
    energies = np.array([*np.ndindex(9, 9)])
    calibration = survey.routes.calibrate(survey.routes.codebook()[0])[0]
    active = survey.routes.graph(board, calibration)
    station_graph = survey.bulk.geometry._atlas([board])[0][0]
    upper = np.triu_indices(board.cells, 1)
    readings = {"native": [], "reversed_second": []}
    held = {}
    rng = np.random.default_rng(28122868)
    for i, j in combinations(range(board.cells), 2):
        times = np.array(
            [
                shortest_path(
                    active(i, j, k + 1, w + 1)[0],
                    unweighted=True,
                    indices=np.arange(0, board.n, 12),
                )[:, ::12][upper]
                for k, w in np.ndindex(9, 9)
            ],
            dtype=np.int16,
        )
        extended = np.column_stack([times, price[energies]])
        direct = np.flatnonzero((upper[0] == i) & (upper[1] == j))
        near = [
            h
            for h, (u, v) in enumerate(zip(*upper, strict=True))
            if (i in (u, v) or j in (u, v)) and station_graph.has_edge(u, v)
        ]
        for tag, stocks in (("native", weights[:3]), ("reversed_second", turned)):
            # Lower address holds the unreversed slot; both orders are
            # physically distinct. Every pair has the same primary writer.
            ref, fresh, shared = stocks
            delta = shared - fresh
            clock, _ = partition(times, ref, fresh, delta)
            both, labels = partition(extended, ref, fresh, delta)
            readings[tag].append(
                {
                    "sites": [i, j],
                    "clock": clock,
                    "clock_and_energy": both,
                    "direct_and_energy": partition(
                        np.column_stack([times[:, direct], price[energies]]),
                        ref,
                        fresh,
                        delta,
                    )[0],
                    "neighbor_times_and_energy": partition(
                        np.column_stack([times[:, near], price[energies]]),
                        ref,
                        fresh,
                        delta,
                    )[0],
                }
            )
            if clock["sum_absolute_count_change"] == 0 and tag not in held:
                groups = np.bincount(labels, weights=delta.ravel()).astype(int)
                mask = groups[labels] > 0
                probabilities = [
                    law.ravel()[mask].sum() / 2592 for law in (fresh, shared)
                ]
                n = 2**26
                recorded = [int(rng.binomial(n, v)) for v in probabilities]
                error = np.sqrt(np.log(2 * 4 / 0.01) / (2 * n))
                diff = (recorded[1] - recorded[0]) / n
                held[tag] = {
                    "sites": [i, j],
                    "starts_per_arm": n,
                    "joint_event_counts": recorded,
                    "positive_joint_classes": np.flatnonzero(groups > 0).tolist(),
                    "joint_event_change_interval": [diff - 2 * error, diff + 2 * error],
                    "change_of_each_joint_class": groups.tolist(),
                    "equivalence_class_by_word": labels.tolist(),
                }
    base, product, _ = weights[:3]

    def entropy(p):
        return float(-np.sum(p / 2592 * np.log(p / 2592)))

    return {
        "full_word_entropy_slope": float(
            -((product - base) / 2592).ravel() @ np.log(base.ravel() / 2592)
        ),
        "full_word_finite_entropy_change": entropy(product) - entropy(base),
        "placements": readings,
        "local_reader": local_reader(
            calibration,
            price,
            (("native", weights[:3]), ("reversed_second", turned)),
            rng,
        ),
        "held": held,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "PAIR_CLOCK_INFORMATION.md",
                "paired_ranging.json",
                "word_routes.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    z = run()
    target = HERE / "pair_clock_information.json"
    if args.write:
        target.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(target.read_bytes()), z)
    for tag, rows in z["placements"].items():
        print(
            tag,
            [
                (k, sum(r[k]["sum_absolute_count_change"] == 0 for r in rows))
                for k in ("clock", "clock_and_energy")
            ],
            z["held"].get(tag, {}).get("sites"),
        )
