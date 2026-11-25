"""Resolve the contact pattern on the same launch as the final source reading."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import beta

from .. import bare_receivers, wall_flight, walls
from ..verify import compare

EXPOSURES = (0, 2, 8, 32)
SHOTS = 2**26
SEED = 81138
ARMS = ((8, False), (8, True), (0, False))
TAIL = 0.01 / (2 * len(EXPOSURES) * len(ARMS) * 2 * 5)
PATH = Path(__file__).with_name("results.json")


def interval(n, total):
    if total == 0:
        return None
    return [
        float(beta.ppf(TAIL, n, total - n + 1)) if n else 0.0,
        float(beta.isf(TAIL, n + 1, total - n)) if n < total else 1.0,
    ]


def contact_masks(anchors):
    bits = np.array(list(product((0, 1), repeat=10)))
    ports = np.stack([bits, 1 - bits], axis=2).reshape(-1, 20)
    edges = [(i, (i + 1) % 8) for i in (1, 3, 5, 7)]
    edges += [(anchors[0], 8), (anchors[1], 14)]
    edges += [(9, 10), (11, 12), (15, 16), (17, 18)]
    mismatch = np.column_stack([ports[:, a] != ports[:, b] for a, b in edges])
    ring_clear = ~mismatch[:, :4].any(axis=1)
    attachments = mismatch[:, 4:6].sum(axis=1)
    left = mismatch[:, 6:8].sum(axis=1)
    right = mismatch[:, 8:10].sum(axis=1)
    primary = ring_clear & (attachments == 0) & (left == 1) & (right == 1)
    distal = primary & mismatch[:, 7] & mismatch[:, 9]
    return mismatch, ring_clear, primary, distal


def initial_states(anchors):
    ref = walls.reference(0)
    inputs = np.asarray(walls.loads(ref))[anchors]
    ground = np.outer(*(walls.wing_ref(8, int(x), False, 3)[2] for x in inputs))
    states = np.zeros((16, 64, 2))
    basis = list(product((0, 1), repeat=4))
    for wall in (0, 1):
        src = tuple(x ^ int(wall and i >= 2) for i, x in enumerate(ref))
        states[basis.index(src), :, wall] = ground.ravel()
    return states.reshape(1024, 2)


def record(state, masks, rng):
    probabilities = np.abs(state) ** 2
    np.testing.assert_allclose(probabilities.sum(), 1, rtol=0, atol=2e-9)
    # Round tiny integration differences before sampling in a fixed word order.
    weights = probabilities.round(11)
    weights /= weights.sum()
    counts = rng.multinomial(SHOTS, weights)
    mismatch, clear, primary, distal = masks
    extinct = int(counts[clear].sum())
    hits, far = (int(counts[m].sum()) for m in (primary, distal))
    return dict(
        counts_by_word=counts.tolist(),
        ordered=extinct,
        primary=hits,
        distal=far,
        intervals=dict(
            ordered=interval(extinct, SHOTS),
            primary=interval(hits, SHOTS),
            distal=interval(far, SHOTS),
            primary_given_ordered=interval(hits, extinct),
            distal_given_ordered=interval(far, extinct),
        ),
        mismatches_on_ordered_words=(counts[clear] @ mismatch[clear]).tolist(),
    )


def no_wall_bound(h, initial, mismatch):
    # Zero mismatches fixes the full load word up to simultaneous reversal.
    zero = ~mismatch.any(axis=1)
    assert zero.sum() == 2
    diagonal = h.diagonal()
    np.testing.assert_allclose(diagonal[zero], 325.0)
    coupling = (h - sparse.diags(diagonal))[zero].toarray()
    np.testing.assert_allclose(coupling @ coupling.T, 2.5 * np.eye(2))
    mean = float(initial @ (h @ initial))
    deviation = np.linalg.norm(h @ initial - mean * initial)
    np.testing.assert_allclose(deviation, 1.0, rtol=0, atol=1e-11)
    assert 62.0 <= mean - 325.0 <= 65.0
    return dict(
        ordered_word_indices=np.flatnonzero(zero).tolist(),
        initial_mean=mean,
        initial_standard_deviation=float(deviation),
        conservative_probability_bound=((1.0 + np.sqrt(2.5)) / 62.0) ** 2,
    )


def measure():
    ring, anchors, wings = wall_flight.extended_geometry()
    initial = initial_states(anchors)
    masks = contact_masks(anchors)
    rng = np.random.default_rng(SEED)
    cases = []
    bound = None
    for k, cut in ARMS:
        h = sparse.csr_matrix(walls.operators(k, cut, 3))
        mask = masks[0].copy()
        if cut:
            mask[:, 4:6] = False
        # Each of the ten dimers has the same diagonal exchange and contrast term.
        np.testing.assert_allclose(
            h.diagonal(), (4 * k + 0.5) * (10 + mask.sum(axis=1))
        )
        if k == 8 and not cut:
            bound = no_wall_bound(h, initial[:, 1], mask)
        flights, _ = bare_receivers.flights(h, initial, EXPOSURES)
        readings = []
        for time, state in zip(EXPOSURES, flights, strict=True):
            readings.append(
                dict(
                    exposure=time,
                    sham=record(state[:, 0], masks, rng),
                    walls=record(state[:, 1], masks, rng),
                )
            )
        cases.append(dict(contrast=k, attachments_cut=cut, readings=readings))
    return dict(
        shots=SHOTS,
        seed=SEED,
        ring=ring,
        anchors=anchors,
        wings=wings,
        word_encoding="big-endian dimer bits; 0 is 02, 1 is 20; ring,left,right",
        all_exposure_wall_launch_bound=bound,
        cases=cases,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = measure()
    if args.write:
        PATH.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n")
    else:
        compare(data, json.loads(PATH.read_text()), "contact_flight")
        print("Verified all joint source and guide histograms.")


if __name__ == "__main__":
    main()
