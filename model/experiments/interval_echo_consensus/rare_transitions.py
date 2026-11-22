"""Independent stratified trials of one physical reversible repair epoch at the exact one-offer reference.

Fixed before acquisition. Every row starts in the named sharp mode with empty scratch banks.
Three fair bits select idle or one of four alternative modes. Every comparison draws fresh
endpoint preparations and actual offer outcomes, keeping its saturation and loss records.
The equilibrium calculation is confined to the analysis following acquisition.
"""

import argparse
from collections import Counter
import json
from math import asin, isclose, log1p, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell

HERE = Path(__file__).resolve().parent

SEED = 57902173
STARTS = (2**26, 2**23, 2**23, 2**18, 2**18)
BATCH = 65536
CAP = 128
MODES = ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1))
OFFER = np.array((0, 0, -1, 1), dtype=np.int8)
RESULT = HERE / "rare_transitions.json"


def source_probability():
    # Symmetric Perron condition for surviving local ternary offers gives b0:b1:b2 = 1:sqrt(2):1.
    middle = sqrt(2) - 1
    angle = 2 * asin(sqrt(middle))
    exchange, _ = law.local_term_parts()
    state = np.eye(law.DIM)[:, law.SRC.index((1, 0))]
    read_port = law.SRC.index((0, 1))
    actual = float(abs((expm(-1j * angle * exchange) @ state)[read_port]) ** 2)
    if abs(actual - middle) > 1e-14:
        raise ValueError("Dimer preparation did not give the required middle outcome.")
    return actual, angle


def endpoints(rng, n, middle_threshold):
    raw = rng.bit_generator.random_raw(2 * n).reshape(-1, 2)
    # Disjoint high 53 bits for the nontrivial dimer read and low bits for its
    # second fair outer rail and the two fair offer bits. Both idle words survive.
    start = np.where((raw >> 11) < middle_threshold, 1, 2 * (raw & 1)).astype(np.int8)
    offers = ((raw >> 1) & 3).astype(np.int8)
    end = start + OFFER[offers]
    valid = np.all((end >= 0) & (end <= 2), axis=1)
    signed = end[:, 0] - end[:, 1]
    # Each code preserves the two starting trits and the two offer words.
    words = (
        (4 * start[:, 0] + offers[:, 0]) * 12 + 4 * start[:, 1] + offers[:, 1]
    ).astype(np.uint8)
    return valid, signed, words


def acquire_row(origin, trials, rng, threshold, witnesses):
    terminal = np.zeros(5, dtype=np.int64)
    proposals = np.zeros(5, dtype=np.int64)
    costs = Counter(
        {
            key: 0
            for key in (
                "idle",
                "accepted",
                "current_wins",
                "cap_expires",
                "pair_draws",
                "pair_saturates",
                "valid_other",
                "match_proposed",
                "match_current",
            )
        }
    )
    examples = {}
    options = [j for j in range(5) if j != origin]
    magnitude, parity = MODES[origin]
    for offset in range(0, trials, BATCH):
        size = min(BATCH, trials - offset)
        selector = (rng.bit_generator.random_raw(size) & 7).astype(np.uint8)
        idle = int(np.count_nonzero(selector < 4))
        terminal[origin] += idle
        costs["idle"] += idle
        for choice, target in enumerate(options, 4):
            accepted = np.flatnonzero(selector == choice)
            proposed_count = accepted.size
            proposals[target] += proposed_count
            new_magnitude, new_parity = MODES[target]
            fields = (
                [(magnitude, new_magnitude)] * witnesses
                if magnitude != new_magnitude
                else []
            ) + ([(parity, new_parity)] if parity != new_parity else [])
            capture_example = (magnitude == 2) != (new_magnitude == 2) and str(
                target
            ) not in examples
            if capture_example:
                attempts_by_gate = np.zeros((size, len(fields)), np.uint8)
                endpoint_words = np.zeros_like(attempts_by_gate)
            for gate, (old, new) in enumerate(fields):
                if accepted.size == 0:
                    break
                pending = accepted
                winners = []
                for attempt in range(1, CAP + 1):
                    if pending.size == 0:
                        break
                    valid, signed_difference, words = endpoints(
                        rng, len(pending), threshold
                    )
                    win = valid & (signed_difference == new)
                    lose = valid & (signed_difference == old)
                    wait = ~(win | lose)
                    wins = int(win.sum())
                    losses = int(lose.sum())
                    valid_pairs = int(valid.sum())
                    costs["pair_draws"] += len(pending)
                    costs["pair_saturates"] += len(pending) - valid_pairs
                    costs["valid_other"] += valid_pairs - wins - losses
                    costs["match_proposed"] += wins
                    costs["match_current"] += losses
                    costs["current_wins"] += losses
                    if wins:
                        winning = pending[win]
                        winners.append(winning)
                        if capture_example:
                            attempts_by_gate[winning, gate] = attempt
                            endpoint_words[winning, gate] = words[win]
                    pending = pending[wait]
                costs["cap_expires"] += len(pending)
                accepted = (
                    np.concatenate(winners) if winners else np.empty(0, dtype=int)
                )
            count = len(accepted)
            terminal[target] += count
            costs["accepted"] += count
            terminal[origin] += proposed_count - count
            if count and capture_example:
                i = int(np.min(accepted))
                examples[str(target)] = dict(
                    epoch=offset + i,
                    initial=origin,
                    proposed=target,
                    compared_signed_words=fields,
                    raw_pairs_through_decision=attempts_by_gate[i].tolist(),
                    winning_endpoint_word_codes=endpoint_words[i].tolist(),
                )
    if (
        terminal.sum() != trials
        or sum(costs[k] for k in ("idle", "accepted", "current_wins", "cap_expires"))
        != trials
    ):
        raise ValueError("Lost a completed or silent epoch.")
    if (
        sum(
            costs[k]
            for k in (
                "pair_saturates",
                "valid_other",
                "match_proposed",
                "match_current",
            )
        )
        != costs["pair_draws"]
    ):
        raise ValueError("Lost a raw pair record.")
    return dict(
        starts=trials,
        proposals=proposals.tolist(),
        terminal=terminal.tolist(),
        counts=dict(costs),
        examples=examples,
    )


def interval(k, n, tail):
    return (
        float(beta.ppf(tail, k, n - k + 1)) if k else 0.0,
        float(beta.ppf(1 - tail, k + 1, n - k)) if k < n else 1.0,
    )


def theoretical(middle, witnesses):
    outer = (1 - middle) / 2
    initial = np.array((outer, middle, outer))
    survivor = (
        (np.diag((2, 2, 2)) + np.diag((1, 1), -1) + np.diag((1, 1), 1)) @ initial
    ) / 4
    s = float(survivor.sum()) ** 2
    q = survivor / survivor.sum()
    p = np.array((q @ q, q[0] * q[1] + q[1] * q[2], q[0] * q[2]))
    target = np.array([p[m] ** witnesses * p[a] for m, a in MODES])
    target /= target.sum()
    P = np.zeros((5, 5))
    cost = np.zeros(5)
    for i, (m, a) in enumerate(MODES):
        for j, (n, b) in enumerate(MODES):
            if i == j:
                continue
            advance = 1 / 8
            for prev, new in ([(m, n)] * witnesses if m != n else []) + (
                [(a, b)] if a != b else []
            ):
                hit = s * (p[prev] + p[new])
                resolve = -np.expm1(CAP * log1p(-hit))
                cost[i] += advance * resolve / hit
                advance *= resolve * p[new] / (p[new] + p[prev])
            P[i, j] = advance
        P[i, i] = 1 - P[i].sum()
    np.testing.assert_allclose(target[:, None] * P, (target[:, None] * P).T, atol=1e-15)
    return s, target, P, cost


def validate_examples(rows, witnesses):
    for origin, r in enumerate(rows):
        for dest, record in r["examples"].items():
            target = int(dest)
            m, a = MODES[origin]
            n, b = MODES[target]
            fields = ([(m, n)] * witnesses if m != n else []) + (
                [(a, b)] if a != b else []
            )
            if tuple(map(tuple, record["compared_signed_words"])) != tuple(fields):
                raise ValueError(
                    "The selected witness coordinates differ from the service."
                )
            for code, attempt, (_, required) in zip(
                record["winning_endpoint_word_codes"],
                record["raw_pairs_through_decision"],
                fields,
                strict=True,
            ):
                left, right = divmod(code, 12)
                x, op_x = divmod(left, 4)
                y, op_y = divmod(right, 4)
                x += int(OFFER[op_x])
                y += int(OFFER[op_y])
                if not (
                    1 <= attempt <= CAP
                    and 0 <= x <= 2
                    and 0 <= y <= 2
                    and x - y == required
                ):
                    raise ValueError(
                        "Archived local histories do not witness the proposed transition."
                    )


def analyze(rows, middle, witnesses):
    counts = np.array([r["terminal"] for r in rows])
    sizes = counts.sum(axis=1)
    measured = counts / sizes[:, None]
    validate_examples(rows, witnesses)
    surv, target, truth, cost = theoretical(middle, witnesses)
    # Simultaneous probability and stationarity-ratio intervals use independent one-epoch
    # preparations in each row. No long-run dwell is counted as an independent discovery.
    tail = 0.05 / (2 * 8)
    ratios = []
    bounds = []
    for i in range(1, 5):
        u = interval(counts[0, i], sizes[0], tail)
        v = interval(counts[i, 0], sizes[i], tail)
        ratios.append(measured[0, i] / measured[i, 0])
        bounds.append((u[0] / v[1], u[1] / v[0]))
    ratios = np.array([1.0, *ratios])
    bounds = np.array([(1.0, 1.0), *bounds])
    direct_pi = ratios / ratios.sum()
    lo = bounds[:, 0] / (bounds[:, 0] + bounds[:, 1].sum() - bounds[:, 1])
    hi = bounds[:, 1] / (bounds[:, 1] + bounds[:, 0].sum() - bounds[:, 0])
    h = measured.T - np.eye(5)
    h[-1] = 1
    empirical_pi = np.linalg.solve(h, np.array([0.0, 0.0, 0.0, 0.0, 1.0]))
    output = dict(
        native_pair_survival=surv,
        stationary_oracle=target.tolist(),
        kernel_oracle=truth.tolist(),
        stationarity_from_empirical_kernel=empirical_pi.tolist(),
        empirical_flow_into_high_if_stationary=float(
            empirical_pi[:3] @ measured[:3, 3:].sum(axis=1)
        ),
        exact_flow_into_high_if_stationary=float(
            target[:3] @ truth[:3, 3:].sum(axis=1)
        ),
        mean_raw_pairs_per_stationary_epoch=float(target @ cost),
        raw_pairs_per_stationary_entry_oracle=float(
            target @ cost / (target[:3] @ truth[:3, 3:].sum(axis=1))
        ),
        direct_probability_estimates=direct_pi.tolist(),
        simultaneous_95_population_bounds=np.column_stack((lo, hi)).tolist(),
        held_even_multiplier=ratios[3] / ratios[1] ** 4,
        held_odd_multiplier=ratios[4] / (ratios[1] ** 3 * ratios[2]),
        held_even_interval=[
            bounds[3, 0] / bounds[1, 1] ** 4,
            bounds[3, 1] / bounds[1, 0] ** 4,
        ],
        held_odd_interval=[
            bounds[4, 0] / (bounds[1, 1] ** 3 * bounds[2, 1]),
            bounds[4, 1] / (bounds[1, 0] ** 3 * bounds[2, 0]),
        ],
        distinct_entries_into_high=int(counts[:3, 3:].sum()),
        exits_from_prepared_high=int(counts[3:, :3].sum()),
    )
    # These are seeded consistency checks, using independently constructed rates after acquisition.
    variances = sizes[:, None] * truth * (1 - truth)
    output["largest_marginal_binomial_standardized_residual"] = float(
        np.max(abs(counts - sizes[:, None] * truth) / np.sqrt(variances))
    )
    return output


def run():
    prob, angle = source_probability()
    threshold = int(prob * 2**53)
    board = cell()
    witness_counts = {
        len({tuple(sorted((p, v))) for p in e for v in board.neighbours[p]} - {e})
        for e in board.seams
    }
    if len(witness_counts) != 1:
        raise ValueError("The source cell has inequivalent seams.")
    (m,) = witness_counts
    seeds = np.random.SeedSequence(SEED).spawn(len(MODES))
    rows = [
        acquire_row(i, n, np.random.Generator(np.random.PCG64(seed)), threshold, m)
        for i, (n, seed) in enumerate(zip(STARTS, seeds, strict=True))
    ]
    return dict(
        seed=SEED,
        starts=list(STARTS),
        batch=BATCH,
        cap=CAP,
        modes=MODES,
        witnesses=m,
        middle_angle=angle,
        middle_probability=prob,
        native_probability_53bit_error=abs(prob - threshold / 2**53),
        acquisition=rows,
        analysis=analyze(rows, prob, m),
    )


def compare(actual, expected, where="receipt"):
    """Require every seed, integer count and material example, with roundoff allowance for analysis."""
    if type(actual) is not type(expected):
        raise ValueError(f"{where} has a different type")
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError(f"{where} has different fields")
        for key in expected:
            compare(actual[key], expected[key], f"{where}.{key}")
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(f"{where} has a different length")
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare(a, b, f"{where}[{i}]")
    elif isinstance(expected, float):
        if not isclose(actual, expected, rel_tol=2e-11, abs_tol=2e-15):
            raise ValueError(f"{where} differs from the archived numerical assessment")
    elif actual != expected:
        raise ValueError(f"{where} differs from the archived record")


def coupling_bootstrap(rows):
    """Exploratory marginal 95% percentile intervals from measured one-epoch transition rows.

    This readout was chosen after acquisition. It uses all rows and complements the
    conservative simultaneous direct-balance intervals already stored in the receipt.
    It neither alters the counts nor supplies a condition for accepting an acquisition.
    """
    counts = np.array([row["terminal"] for row in rows])
    starts = counts.sum(axis=1)
    measured = counts / starts[:, None]
    rng = np.random.default_rng(3961405)
    draws = [
        rng.multinomial(int(n), p, size=8192) / n
        for n, p in zip(starts, measured, strict=True)
    ]
    kernels = np.concatenate((measured[None, :, :], np.stack(draws, axis=1)))
    equations = kernels.transpose(0, 2, 1) - np.eye(5)
    equations[:, -1, :] = 1.0
    rhs = np.zeros((len(kernels), 5, 1))
    rhs[:, -1, 0] = 1.0
    pi = np.linalg.solve(equations, rhs)[:, :, 0]
    beta_read = np.log(pi[:, 1] / pi[:, 2])
    samples = {
        "beta": beta_read,
        "magnitude_one_energy": np.log(pi[:, 0] / pi[:, 1]) / beta_read,
        "magnitude_two_even_energy": np.log(pi[:, 0] / pi[:, 3]) / beta_read,
    }
    return {
        key: [float(value[0]), *map(float, np.quantile(value[1:], [0.025, 0.975]))]
        for key, value in samples.items()
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", action="store_true")
    args = cli.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.record:
        with RESULT.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text()))
    a = result["analysis"]
    print(
        a["distinct_entries_into_high"],
        a["exits_from_prepared_high"],
        a["held_even_multiplier"],
        a["held_odd_multiplier"],
    )
    print(
        "Exploratory row-bootstrap estimates and marginal 95% intervals:",
        coupling_bootstrap(result["acquisition"]),
    )
