"""Independent first-destination counts from protected link-phase tickets."""

import math

import numpy as np
from scipy.sparse import csr_matrix


def acquire(step, start, pins, seed, copies, *, excursion=False):
    n = len(step)
    a, b = pins
    live = step.copy()
    live[[a, b]] = 0
    live[a, a] = live[b, b] = 12
    # Integer counts of all native slot and phase paths in eight appointments.
    # Draw a multinomial per occupied start, exactly as for independent trials.
    propagator = csr_matrix(live)
    block = np.eye(n, dtype=np.int64)
    for _ in range(8):
        block = propagator @ block
    den = 12**8
    assert np.all(block.sum(1) == den)
    destinations = [np.flatnonzero(row) for row in block]
    probabilities = [row[z] / den for row, z in zip(block, destinations, strict=True)]
    counts = np.zeros(n, dtype=np.int64)
    rng = np.random.default_rng(seed)
    if excursion:
        counts = rng.multinomial(copies, step[start] / 12)
    else:
        counts[start] = copies
    arrived = counts[[a, b]].copy()
    counts[[a, b]] = 0
    for k in range(1, 12001):
        nxt = np.zeros_like(counts)
        for i in np.flatnonzero(counts):
            z = destinations[i]
            nxt[z] += rng.multinomial(counts[i], probabilities[i])
        arrived += nxt[[a, b]]
        nxt[[a, b]] = 0
        counts = nxt
        if counts.sum() == 0:
            break
    assert arrived.sum() + counts.sum() == copies
    return {
        "first": arrived.tolist(),
        "pending": int(counts.sum()),
        "trials": copies,
        "blocks": k,
    }


def source_and_receiver(
    steps,
    probabilities,
    exact,
    observed,
    initial_contrast,
    copies_per_start,
    launch,
    pins,
):
    copies = 2**26
    rows = [acquire(s, launch, pins, 418973 + j, copies) for j, s in enumerate(steps)]
    p = np.array([z["first"][1] / copies for z in rows])
    pending = np.array([z["pending"] / copies for z in rows])
    obs = np.array(observed)
    # Conditional on the separately drawn stratified source bank the remote
    # observations have known weights. Its finite-count width uses the full
    # spread of the conditional receiver probabilities, bounded simultaneously
    # from the remote counts. The launch was selected using calculated curves.
    spread = (
        np.ptp(p)
        + max(pending)
        + 2 * math.sqrt(math.log(2 * len(rows) / 0.002) / (2 * copies))
    )
    remote = (
        math.sqrt((obs @ obs) * math.log(2 / 0.003) / (2 * copies)) + abs(obs) @ pending
    )
    local = spread * math.sqrt(
        (initial_contrast @ initial_contrast)
        * math.log(2 / 0.003)
        / (2 * copies_per_start)
    )
    estimate = float(obs @ (p - p[0]))
    # Each start stratum has the same exact number of path records, hence the
    # signed source estimate has zero total. Float roundoff is included here.
    pad = remote + local + 1e-10
    return {
        "transfers": rows,
        "destination_probabilities": p.tolist(),
        "receiver_width": remote,
        "source_width": local,
        "change": estimate,
        "simultaneous_interval": [estimate - pad, estimate + pad],
        "coverage_at_least": 0.992,
        "native_ticket_prediction": float(exact @ (probabilities - probabilities[0])),
    }
