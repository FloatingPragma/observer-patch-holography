"""Recover a common angular quadrature from short, labelled material flights."""

import json
from fractions import Fraction
from math import sin

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.stats import beta

from model.candidate_laws.interval_echo.clock import _dimer_pulse

from .. import measure as tracer

DWELLS = (Fraction(1, 512), Fraction(1, 256))
SHOTS = 2**34


def exchange_certificate(twice, tag):
    """Uniform rotating-frame bound for the entire resident subspace at each tag launch."""
    hopping = csr_matrix(np.where(tag[:, None] != tag[None, :], twice, 0))
    internal = csr_matrix(twice) - hopping
    commutator = internal @ hopping - hopping @ internal
    bound = Fraction(
        int(abs(commutator).sum(axis=1).max())
        + int(abs(hopping).sum(axis=1).max()) ** 2,
        4,
    )
    receivers = np.unique(tag)
    edges = np.zeros((len(receivers), len(receivers)), dtype=bool)
    for a in receivers:
        for b in receivers:
            block = hopping[tag == b][:, tag == a]
            if block.nnz:
                assert np.all(block.data == -1)
                np.testing.assert_array_equal(
                    abs(block).sum(axis=0), np.ones((1, block.shape[1]))
                )
                np.testing.assert_array_equal(
                    abs(block).sum(axis=1), np.ones((block.shape[0], 1))
                )
                edges[a, b] = True
    return bound, edges


def bearing_frame(edges):
    """The slow projector and the six bearing dyads from an anonymous adjacency table."""
    assert np.array_equal(edges, edges.T)
    steps = shortest_path(csr_matrix(edges), directed=False, unweighted=True)
    near = (steps == 1).astype(int)
    far = (steps == 2).astype(int)
    opposite = (steps == 3).astype(int)
    np.testing.assert_array_equal(opposite.sum(axis=0), np.ones(12))
    even = np.eye(12, dtype=int) - opposite
    odd = near - far
    np.testing.assert_array_equal(5 * even @ even + odd @ odd, 20 * even)
    np.testing.assert_array_equal(even @ odd + odd @ even, 4 * odd)
    np.testing.assert_array_equal(near @ even, odd)
    np.testing.assert_array_equal(near @ odd, 5 * even)
    assert np.trace(even) == 12 and np.trace(odd) == 0
    gram = even + odd / np.sqrt(5)
    axes = [
        (a, int(np.argmax(opposite[a])))
        for a in range(12)
        if a < np.argmax(opposite[a])
    ]
    positive = [a for a, _ in axes]
    dyads = gram[np.ix_(positive, positive)] ** 2
    np.testing.assert_allclose(5 * dyads, np.ones((6, 6)) + 4 * np.eye(6), atol=1e-13)
    return axes, gram


def simulate():
    twice, basis, rim, outside = tracer.source(8)
    tag = np.array([v for v, _ in basis])
    pairs = [w for _, w in basis]
    eig, vecs = eigh(twice / 2)
    bound, expected = exchange_certificate(twice, tag)
    ref = sorted({w for w in pairs if outside not in w})
    rng = np.random.default_rng(61024)
    labels = rng.permutation(12)
    exposures = []
    reference_rng = np.random.default_rng(40957)
    settings = 2 * (3 * 12 * 12 + 1)
    tail = 0.01 / (2 * settings)
    for dt in DWELLS:
        remainder = dt**2 * bound / 2
        threshold = dt**2 / 16
        assert (
            remainder**2 < threshold < (dt / 2 - remainder) ** 2 and dt / 2 > remainder
        )
        # An isolated material exchange at the same shutter gives an independent
        # optical clock standard on the registered short branch.
        hit_rate = abs(_dimer_pulse(float(dt))[1, 0]) ** 2
        standard_hits = int(reference_rng.binomial(SHOTS, hit_rate))
        measured_lo = beta.ppf(tail, standard_hits, SHOTS - standard_hits + 1)
        measured_hi = beta.ppf(1 - tail, standard_hits + 1, SHOTS - standard_hits)
        time_lo, time_hi = 2 * np.arcsin(np.sqrt([measured_lo, measured_hi]))
        worst_remainder = float(bound) * time_hi**2 / 2
        assert (
            worst_remainder**2
            < measured_lo / 4
            < measured_hi / 4
            < (time_lo / 2 - worst_remainder) ** 2
        )
        assert time_lo > 2 * worst_remainder
        wave = (vecs * np.exp(-1j * float(dt) * eig)) @ vecs.T
        tag_probs = np.array(
            [np.sum(abs(wave[tag == a]) ** 2, axis=0) for a in range(12)]
        )
        readouts = {}
        for name, stock in (
            ("uniform_ball", ref),
            ("split", tracer.clock_ball.SPLIT),
            ("spread", tracer.clock_ball.SPREAD),
        ):
            selected = np.array([w in stock for w in pairs])
            banks = []
            for a in labels:
                joint = tag_probs[:, selected & (tag == a)].sum(axis=1) / len(stock)
                admit = sum(a not in w for w in stock) / len(stock)
                np.testing.assert_allclose(joint.sum(), admit, atol=1e-12)
                # Preserve blocked launches. Normalize the admitted conditional vector
                # before asking the multinomial sampler for its disjoint terminal bins.
                law = np.r_[joint[labels] / joint.sum() * admit, 1 - admit]
                banks.append(rng.multinomial(SHOTS, law))
            banks = np.array(banks)
            hits = banks[:, :12]
            admitted = hits.sum(axis=1)[:, None]
            lo = beta.ppf(tail, np.maximum(hits, 1), admitted - hits + 1)
            lo[hits == 0] = 0.0
            hi = beta.ppf(1 - tail, hits + 1, np.maximum(admitted - hits, 1))
            hi[hits == admitted] = 1.0
            edges = lo > measured_hi / 4
            np.fill_diagonal(edges, False)
            unclassified = (lo <= measured_hi / 4) & (hi >= measured_lo / 4)
            np.fill_diagonal(unclassified, False)
            assert not unclassified.any()
            axes, gram = bearing_frame(edges)
            np.testing.assert_array_equal(edges, expected[np.ix_(labels, labels)])
            nonedge = ~edges & ~np.eye(12, dtype=bool)
            readouts[name] = dict(
                counts=banks.tolist(),
                recovered_edges=int(edges.sum() // 2),
                largest_nonedge_upper=float(hi[nonedge].max()),
                smallest_edge_lower=float(lo[edges].min()),
                antipodes_in_receiver_labels=axes,
                gram_rank=int(np.linalg.matrix_rank(gram)),
            )
        exposures.append(
            dict(
                dwell=str(dt),
                standard_exchange_hits=standard_hits,
                clock_interval=[float(time_lo), float(time_hi)],
                empirical_threshold_interval=[
                    float(measured_lo / 4),
                    float(measured_hi / 4),
                ],
                decision_probability=float(threshold),
                universal_nonedge_upper=str(remainder**2),
                universal_edge_lower=str((dt / 2 - remainder) ** 2),
                counts=readouts,
            )
        )
    control, _, _, _ = tracer.source(0)
    comparison_bound, control_edges = exchange_certificate(control, tag)
    np.testing.assert_array_equal(expected, control_edges)
    return dict(
        word_count=len(basis),
        shots_per_launch=SHOTS,
        anonymous_labels=labels.tolist(),
        rotating_frame_bound=str(bound),
        zero_contrast_bound=str(comparison_bound),
        zero_exchange_offsite_probability=0,
        isolated_contact_at_first_dwell=sin(float(DWELLS[0]) / 2) ** 2,
        six_dyad_spectrum=["4/5"] * 5 + ["2"],
        exposures=exposures,
    )


if __name__ == "__main__":
    print(json.dumps(simulate(), indent=2, allow_nan=False))
