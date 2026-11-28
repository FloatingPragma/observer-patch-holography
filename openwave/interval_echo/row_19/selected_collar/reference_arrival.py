"""Transmit the measured opposite-reference modular score by comparing ordinary records."""

from fractions import Fraction as F
from math import fsum, log, sqrt

from mpmath import iv
import numpy as np

from . import coherent_balls as coherent

local = coherent.local
bounds = coherent.balls.bounds
SEED = 28122921
LIMIT = 12288
PER_WORD = 4096


def programs(reference, fresh, shared):
    # All three members of each marginal class are retained separately.
    normal, high, low = (np.flatnonzero(reference[0] == n) for n in (32, 48, 16))

    def cycle(left, right):
        u = np.zeros(9, dtype=int)
        u[left] = 1
        u[right] = -1
        return np.outer(u, u)

    quiet = cycle(normal, high)
    levels = {
        "fresh": fresh,
        "shared": shared,
        "same_serial": fresh + 16 * quiet,
        "opposed_serial": fresh + 32 * quiet - 2 * cycle(normal, low),
    }
    for m in levels.values():
        assert np.min(m) >= 0
        for axis in (0, 1):
            np.testing.assert_array_equal(m.sum(axis=axis), fresh.sum(axis=axis))
    return levels


def panel_lengths(trials, rng):
    # A proposal in bucket j chooses each k there with probability
    # 1/(num_buckets * 2^j). Exactly 2^j of k possible ranks accept.
    lengths = np.empty(trials, dtype=int)
    waiting = np.arange(trials)
    while len(waiting):
        lower = 2 ** rng.integers(LIMIT.bit_length(), size=len(waiting))
        k = lower + rng.integers(lower)
        selected = (k <= LIMIT) & (rng.integers(k) < lower)
        lengths[waiting[selected]] = k[selected]
        waiting = waiting[~selected]
    return lengths


def compare_sharp(symbols, ref, rng):
    # Literal independent reference observations until a match or deadline.
    target = np.repeat(np.arange(len(symbols)), PER_WORD)
    budget = panel_lengths(len(target), rng)
    live = np.arange(len(target))
    success = np.zeros(len(symbols), dtype=int)
    comparisons = 0
    while len(live):
        other = rng.choice(len(symbols), size=len(live), p=ref)
        comparisons += len(live)
        unequal = np.any(symbols[target[live]] != symbols[other], axis=1)
        live = live[unequal]
        budget[live] -= 1
        end = budget[live] == 0
        success += np.bincount(target[live[end]], minlength=len(symbols))
        live = live[~end]
    return success.tolist(), comparisons


def count_record_minima(symbols, ref, rng):
    # Successive insertion ranks are the online coordinates of a uniform
    # random permutation. A zero rank records a new minimum.
    word = np.repeat(np.arange(len(symbols)), PER_WORD)
    live = np.arange(len(word))
    totals = np.zeros(len(symbols), dtype=int)
    comparisons = 0
    for k in range(1, LIMIT + 1):
        if not len(live):
            break
        reference = rng.choice(len(symbols), size=len(live), p=ref)
        comparisons += len(live)
        unequal = np.any(symbols[word[live]] != symbols[reference], axis=1)
        live = live[unequal]
        record = rng.integers(k, size=len(live)) == 0
        totals += np.bincount(word[live[record]], minlength=len(symbols))
    return totals.tolist(), comparisons


def read_settings(weights, prob, serial, rng):
    p = weights.ravel() / 2592
    # Complete native symbols and the two proposed terminal admissions are
    # acquired in independent reset cohorts. Joint distributions of the
    # fair length panel are integrated for the batched large-sample draws.
    counts = coherent.draws(p, rng)
    admission = float(p @ prob)
    return {
        "word_counts": counts,
        "collision_accept_counts": coherent.draws([admission, 1 - admission], rng),
        "serial_accept_counts": coherent.draws([p @ serial, 1 - p @ serial], rng),
        "collision_acceptance": admission,
        "serial_acceptance": str(F(int(weights.ravel() @ serial), 2592)),
    }


def run():
    iv.dps = 50
    rng = np.random.default_rng(SEED)
    ref, fresh, shared = local.survey.standards(local.survey.bulk.federation(0))[0][:3]
    stocks = programs(ref, fresh, shared)
    read = local.inputs("pair_clock_information")["local_reader"]
    letters = read["calibrations"]["fixed_contact"][
        "sharp_records_in_calibration_order"
    ]
    symbols = np.array([a + b for a in letters for b in letters])
    assert len(np.unique(symbols, axis=0)) == len(symbols) == 81
    serial = np.array(
        local.inputs("shared_contact")["serial_route"]["sharp_arrivals"]
    ).ravel()
    assert np.array_equal(serial, (symbols[:, 1] * symbols[:, 3]))
    sizes = np.arange(1, LIMIT + 1)
    scale = sum(1 / iv.mpf(k) for k in range(1, LIMIT + 1))
    refprob = ref.ravel() / 2592
    values = {
        n: fsum(np.exp(sizes * np.log1p(-n / 2592)) / sizes) / float(scale.mid)
        for n in np.unique(ref)
    }
    prob = np.array([values[n] for n in ref.ravel()])
    counts, comparisons = compare_sharp(symbols, refprob, rng)
    minima, comparisons_to_recurrence = count_record_minima(symbols, refprob, rng)
    ref_counts = coherent.draws(refprob, rng)
    half = sqrt(log(2 * 10000 / 0.01) / (2 * coherent.N)) + 1e-13

    def bands(counts):
        return [
            iv.mpf([max(0, k / coherent.N - half), min(1, k / coherent.N + half)])
            for k in counts
        ]

    reference = bands(ref_counts)
    pmin = min(float(z.a) for z in reference) - 1e-16
    p = iv.mpf(pmin)
    # The harmonic normalization times the finite acceptance underestimates each modular cost.
    # All individual errors lie between zero and one common absolute cap;
    # a difference of two normalized mixtures inherits that same cap.
    tail = (1 - p) ** (LIMIT + 1) / ((LIMIT + 1) * p)
    for n, v in values.items():
        true = -iv.log(iv.mpf(int(n)) / 2592)
        assert (
            float(true.a) - float(tail.b) - 2e-14
            <= float(scale.a) * v
            <= float(true.b) + 2e-14
        )
    arms = {k: read_settings(z, prob, serial, rng) for k, z in stocks.items()}
    f = bands(arms["fresh"]["word_counts"])
    base = bands(arms["fresh"]["collision_accept_counts"])[0]
    for name, row in arms.items():
        measured = bands(row["word_counts"])
        q = -sum(
            (a - b) * iv.log(81 * p)
            for a, b, p in zip(measured, f, reference, strict=True)
        )
        change = bands(row["collision_accept_counts"])[0] - base
        exact = -sum(
            iv.mpf(int(v)) / 2592 * iv.log(iv.mpf(int(n)) / 2592)
            for v, n in zip((stocks[name] - fresh).ravel(), ref.ravel(), strict=True)
        )
        row["source_program"] = stocks[name].tolist()
        row["measured_modular_contrast"] = bounds(q)
        row["sharp_law_modular_contrast"] = bounds(exact)
        row["measured_arrival_change"] = bounds(change)
        row["arrival_from_word_measurements"] = bounds(
            (q + iv.mpf([-float(tail.b), float(tail.b)])) / scale
        )
    return {
        "seed": SEED,
        "panel_limit": LIMIT,
        "harmonic_normalization": bounds(scale),
        "record_minima_by_word": minima,
        "reference_comparisons_to_recurrence": comparisons_to_recurrence,
        "reference_word_counts": ref_counts,
        "measured_reference_floor": pmin,
        "max_charge_error_nats": float(tail.b) + 1e-22,
        "literal_starts_per_word": PER_WORD,
        "literal_reference_comparisons": comparisons,
        "literal_passes_by_word": counts,
        "pass_rate_by_word": prob.tolist(),
        "batched_starts": coherent.N,
        "arms": arms,
        "inputs": coherent.steady.exp.save_inputs(
            "ARCHIVED_OVERLAP.md", "pair_clock_information.json", "shared_contact.json"
        ),
    }


if __name__ == "__main__":
    result = run()
    coherent.steady.exp.emit("reference_arrival", result)
    print(result["max_charge_error_nats"])
    for k, v in result["arms"].items():
        print(
            k,
            v["measured_modular_contrast"],
            v["measured_arrival_change"],
            v["serial_acceptance"],
        )
