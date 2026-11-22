"""Draw 20 million finite pair histories in small batches and return the aggregate event ledger.

Run from the repository root with python -m model.experiments.interval_echo_consensus.prospective_offer_capture.
The pseudorandom engine samples the verified native sharp Born transition probabilities.
"""

import hashlib
import json
from itertools import combinations
from math import sqrt
import numpy as np
from scipy.linalg import expm
from scipy.stats import norm

from model.candidate_laws.interval_echo.candidate import native_law as law
from model.executable_model import federation, rotations
from model.experiments.interval_echo_consensus.offer_service import (
    _source_chains,
    _transport,
    signed_probabilities,
)


FROZEN_PROTOCOL = """# A fresh acquisition at seventeen offers

Prospective local capture, 22 November 2026. The 17-offer length is fixed from earlier work. Its previous phase-tolerance decision consulted higher contrast. The following primary forecast tests new simulated draws under a different initial distribution and observer arrangement. It gives no independent validation of that earlier length choice. Exact preliminary source probabilities were available when selecting achievable counts.

Before drawing samples, fix the initial port probability as equal mass at 0 and 2, with no load 1. Every window uses two independent sharp outcomes of half contacts on single units, one fair bit for idle or offer and a second fair bit for direction. Track all idles and the stage and sign of every failed offer for both chains of a pair. The two terminal loads of each surviving pair are read sharply. All rejected and incompletely buffered trials count.

Use NumPy PCG64, chunk size 16,384, 4,000,000 training pair attempts with offer seed 17012026 and parity seed 17012027. Test 16,000,000 new pair attempts with seeds 17022026 and 17022027. Read the training pairs individually, including their native even/odd parity outcomes and a fixed positive chart for labeling signed words. Buffer accepted test pairs in order, four at a time. Three records come from the rotated face (38,39), (38,40), (39,40) of cell 3 on the 20-cell federation. The independent parity reporter uses (36,37) and the negative chart. A single parity observation on a sharp nonzero last-witness word has equal probabilities; a zero word is even. On disagreement the hand is left alone. Compare the chart word only after this parity observation, retaining rejection and public quiet outcomes.

The primary withheld observable is the fraction of all completely buffered four-source cohorts with an accepted |d|=1 hand, pooling its two parities. From signed training frequencies q0, q+ and q−, predict (q+³ + q−³)(q0 + q−)/2. The two orientations are estimated separately before changing chart. Declare agreement within three combined standard errors using the multinomial delta variance for training and the binomial variance for held cohorts. Report the held 95% Wilson interval, training and held standard errors, as well as direct stage counts. All other endpoints are descriptive. In particular the large aperture eight-witness unanimous magnitude-two frequency remains an exact model extrapolation unless directly observed in its own enrollment. Endpoint seeds, sample sizes, geometry, cutoff and test direction are fixed here and shall not be changed after seeing the data.
"""


def layout():
    c = federation(0)
    rotated = min(
        tuple(sorted(r[v] for v in (0, 1, 5)))
        for r in rotations()
        if not {r[v] for v in (0, 1, 5)}.intersection((0, 1, 5))
    )
    verts = tuple(v + 36 for v in rotated)
    stations = [
        *combinations(verts, 2),
        next(e for e in c.seams if all(v // 12 == 3 and v not in verts for v in e)),
    ]
    chains = _source_chains(c, 8, 17, {3})
    used = {tuple(sorted((a, b))) for chain in chains for a, b in zip(chain, chain[1:])}
    pending = {chain[-1] for chain in chains}
    routes = []
    for i, edge in enumerate(stations):
        for k, port in enumerate(edge):
            start = chains[2 * i + k][-1]
            pending.remove(start)
            routes.append(_transport(c, start, port, used, pending))
        used.add(edge)
    return dict(cells=c.cells, stations=stations, source_chains=chains, routes=routes)


def draw_pair_records(attempts, seed):
    draw = np.random.Generator(np.random.PCG64(seed))
    # Each row contains source population before a window, legitimate holds,
    # increases, decreases, saturated low, saturated high and already killed.
    steps = np.zeros((17, 7), dtype=np.int64)
    survivors = []
    branches = np.zeros(3, dtype=np.int64)
    for base in range(0, attempts, 16384):
        n = min(16384, attempts - base)
        x = 2 * draw.integers(2, size=2 * n, dtype=np.int8)
        for stage in range(17):
            live = np.flatnonzero(x >= 0)
            # Values 0,1 encode the two idle words; 2=downward and 3=upward.
            bits = draw.integers(4, size=len(live), dtype=np.int8)
            before = x[live]
            low = (before == 0) & (bits == 2)
            high = (before == 2) & (bits == 3)
            counts = np.bincount(bits, minlength=4)
            bad0, bad2 = np.count_nonzero(low), np.count_nonzero(high)
            steps[stage] += [
                len(live),
                counts[0] + counts[1],
                counts[3] - bad2,
                counts[2] - bad0,
                bad0,
                bad2,
                2 * n - len(live),
            ]
            after = before + (bits == 3).astype(np.int8) - (bits == 2).astype(np.int8)
            after[low | high] = -2
            x[live] = after
        pair = x.reshape(-1, 2)
        success = (pair >= 0).sum(axis=1)
        branches += np.bincount(success, minlength=3)
        found = pair[success == 2]
        survivors.append(found[:, 0] - found[:, 1])
    signs = np.concatenate(survivors)
    if int(branches.sum()) != attempts or len(signs) != branches[2]:
        raise ValueError("Incomplete signed-offer outcomes")
    if not np.all(steps[:, 0] == steps[:, 1:6].sum(axis=1)):
        raise ValueError("Incomplete per-window accounting")
    return signs, dict(
        pair_attempts=attempts,
        zero_one_two_chains_valid=branches.tolist(),
        window_fields=[
            "active",
            "idle",
            "up_valid",
            "down_valid",
            "fail_low",
            "fail_high",
            "previously_failed",
        ],
        windows=steps.tolist(),
        signed_counts_order_minus2_to_plus2=np.bincount(
            signs + 2, minlength=5
        ).tolist(),
    )


def native_fair_controls():
    J = np.eye(law.DIM)[:, : law.N_SRC]
    h0, h1 = law.local_term_parts()
    for k in (0, 1.3, 8):
        half = expm(-1j * np.pi / 2 * (h0 + k * h1))
        initial = J[:, law.SRC.index((1, 0))]
        np.testing.assert_allclose(
            abs(half @ initial) ** 2,
            (initial + J[:, law.SRC.index((0, 1))]) / 2,
            atol=1e-12,
        )
        echo = law.echo()
        phases = np.repeat([1, *(1j if d < 0 else 1 for _, d in law.LABELS)], law.N_SRC)
        analyzer = half @ echo @ (phases[:, None] * (echo @ J))
        mask = np.tile([a <= b for a, b in law.SRC], 1 + law.N_LAB)
        even = (abs(analyzer[mask]) ** 2).sum(axis=0)
        np.testing.assert_allclose(
            even, [1 if a == b else 0.5 for a, b in law.SRC], atol=1e-12
        )


def parity_read(signs, seed):
    bits = np.random.Generator(np.random.PCG64(seed)).integers(
        2, size=len(signs), dtype=np.int8
    )
    return bits * (signs != 0)


def measure():
    mounted = layout()
    native_fair_controls()
    train, a = draw_pair_records(4_000_000, 17012026)
    tb = parity_read(train, 17012027)
    a["native_word_parity_counts"] = (
        np.bincount(2 * (train + 2) + tb, minlength=10).reshape(5, 2).tolist()
    )
    q = np.bincount(train + 2, minlength=5) / len(train)
    minus, zero, plus = q[1:4]
    factor1 = (plus**3 + minus**3) / 2
    factor2 = zero + minus
    pred = factor1 * factor2
    gradient = np.array(
        [0, 1.5 * minus**2 * factor2 + factor1, factor1, 1.5 * plus**2 * factor2, 0]
    )
    se_train = sqrt(gradient @ (np.diag(q) - np.outer(q, q)) @ gradient / len(train))
    hold, b = draw_pair_records(16_000_000, 17022026)
    cohorts = len(hold) // 4
    b["unbuffered_tail"] = hold[4 * cohorts :].tolist()
    four = hold[: 4 * cohorts].reshape(-1, 4)
    equal = (four[:, :3] == four[:, 0, None]).all(axis=1)
    hand = four[equal, 2]
    parity = parity_read(hand, 17022027)
    extra = four[equal, 3]
    accept = extra == -parity
    b["buffered_cohorts"] = cohorts
    b["disagreed_before_parity"] = int(np.count_nonzero(~equal))
    b["agreement_magnitude_parity_counts"] = (
        np.bincount(2 * abs(hand) + parity, minlength=6).reshape(3, 2).tolist()
    )
    b["parity_reject_counts"] = (
        np.bincount(2 * abs(hand[~accept]) + parity[~accept], minlength=6)
        .reshape(3, 2)
        .tolist()
    )
    b["accepted_counts"] = (
        np.bincount(2 * abs(hand[accept]) + parity[accept], minlength=6)
        .reshape(3, 2)
        .tolist()
    )
    y = int(np.count_nonzero(accept & (abs(hand) == 1)))
    obs = y / cohorts
    se_hold = sqrt(obs * (1 - obs) / cohorts)
    z95 = norm.ppf(0.975)
    w = z95 * z95 / cohorts
    mid = (obs + w / 2) / (1 + w)
    rad = z95 * sqrt(obs * (1 - obs) / cohorts + w / (4 * cohorts)) / (1 + w)
    zscore = (obs - pred) / sqrt(se_train**2 + se_hold**2)
    p0, p1, p2 = signed_probabilities(17, (1, 0, 1))
    return dict(
        protocol_sha256=hashlib.sha256(FROZEN_PROTOCOL.encode()).hexdigest(),
        material_layout=mounted,
        training=a,
        held=b,
        enrolled_test=dict(
            forecast=pred,
            training_standard_error=se_train,
            held_events=y,
            held_fraction=obs,
            held_standard_error=se_hold,
            held_wilson95=[mid - rad, mid + rad],
            discrepancy_standard_errors=zscore,
            within_3_standard_errors=bool(abs(zscore) <= 3),
            separate_exact_dp=float(p1**3 * (p0 + p1)),
        ),
        independent_dp_signed=[float(p2), float(p1), float(p0), float(p1), float(p2)],
    )


if __name__ == "__main__":
    print(json.dumps(measure(), indent=2, allow_nan=False))
