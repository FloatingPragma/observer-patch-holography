#!/usr/bin/env python3
"""Tests for the event-manifold reconstruction receipts (GitHub #525)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from event_manifold_reconstruction import (  # noqa: E402
    LocatedRecord,
    bilipschitz_distortion,
    boost_matrix,
    branching_records,
    causal_labels,
    class_positions,
    embed_h3,
    eta_dot,
    event_classes,
    fit_cone_quadric,
    frame_reconstruct,
    generate_records,
    generic_events,
    germs_from_records,
    icosahedral_cap_normals,
    label_inflated_records,
    pair_verdict,
    pairwise_intervals,
    pca_dimension,
    response_matrix,
    signature_of,
    worldline_events,
)

STAGES = 6


def true_chart(events: np.ndarray, taus: np.ndarray) -> np.ndarray:
    return np.stack([np.concatenate(([t], x[1:])) for x, t in zip(events, taus)])


# ---------------------------------------------------------------------------
# frame reconstruction (rank-4 receipt E3)
# ---------------------------------------------------------------------------

def test_frame_is_rank_four_and_reconstructs_exactly():
    normals = icosahedral_cap_normals()
    b = response_matrix(normals)
    assert np.linalg.matrix_rank(b) == 4
    x = embed_h3(np.array([0.3, -0.5, 0.2]))
    h = np.array([eta_dot(x, n) for n in normals])
    x_rec = frame_reconstruct(h, b)
    # h = B x with B rows (-n^0, n^1, n^2, n^3): the left inverse is exact
    assert np.allclose(x_rec, x, atol=1e-9)


# ---------------------------------------------------------------------------
# events, separation, coincidence (E1-E2, thm:event-4-manifold)
# ---------------------------------------------------------------------------

def test_distinct_events_are_certified_separated():
    events, taus = generic_events(6)
    records = generate_records(events, taus, STAGES)
    classes, ambiguous = event_classes(records)
    assert len(classes) == 6
    assert not ambiguous


def test_coincident_germs_merge_to_one_event():
    events, taus = generic_events(5)
    r1 = generate_records(events, taus, STAGES, seed=5)
    r2 = generate_records(events, taus, STAGES, seed=6)
    n = len(events)
    for rec in r2:
        rec.token += n
    classes, ambiguous = event_classes(r1 + r2)
    assert len(classes) == n
    assert not ambiguous
    assert all(len(cls) == 2 for cls in classes)


def test_hausdorff_separation_gap_shrinks_with_refinement():
    events, taus = generic_events(4)
    records = generate_records(events, taus, STAGES)
    germs = germs_from_records(records)
    radii = [germs[0][s].radius for s in range(STAGES)]
    assert all(radii[i] > radii[i + 1] for i in range(STAGES - 1))


# ---------------------------------------------------------------------------
# chart quality and dimension (thm:event-4-manifold)
# ---------------------------------------------------------------------------

def test_chart_is_bilipschitz_and_dimension_four():
    events, taus = generic_events(24)
    records = generate_records(events, taus, STAGES)
    classes, _ = event_classes(records)
    assert len(classes) == 24
    order = [cls[0] for cls in classes]
    rec_pts = class_positions(records, classes)
    tru_pts = true_chart(events, taus)[order]
    assert bilipschitz_distortion(tru_pts, rec_pts) < 1.25
    assert pca_dimension(rec_pts) == 4


# ---------------------------------------------------------------------------
# causal order and signature (thm:event-causal-structure)
# ---------------------------------------------------------------------------

def test_signature_recovered_from_ancestry_labels():
    events, taus = generic_events(40, seed=21)
    records = generate_records(events, taus, STAGES, seed=22)
    classes, _ = event_classes(records)
    order = [cls[0] for cls in classes]
    rec_pts = class_positions(records, classes)
    tru_pts = true_chart(events, taus)[order]
    labels = causal_labels(tru_pts)          # intrinsic ancestry labels
    assert any(flag for *_ij, flag in labels)
    q = fit_cone_quadric(rec_pts, labels)    # fit uses reconstructed chart
    neg, pos = signature_of(q, tol=1e-3 * np.linalg.norm(q, 2))
    assert (neg, pos) == (1, 3)
    # exactly one time direction: negative eigenvector is timelike
    evals, evecs = np.linalg.eigh(q)
    v_time = evecs[:, 0]
    assert eta_dot(v_time, v_time) < 0.0


def test_time_orientation_from_transaction_direction():
    events, taus = generic_events(30, seed=31)
    tru_pts = true_chart(events, taus)
    labels = causal_labels(tru_pts)
    for i, j, causal in labels:
        if causal:
            # ancestry always points forward in the produced clock
            assert tru_pts[j][0] > tru_pts[i][0]


# ---------------------------------------------------------------------------
# invariance (lem:event-quotient-invariance)
# ---------------------------------------------------------------------------

def test_lorentz_boost_gauge_and_schedule_invariance():
    events, taus = generic_events(8, seed=41)
    records = generate_records(events, taus, STAGES, seed=42)
    classes, _ = event_classes(records)
    rec_pts = class_positions(records, classes)
    intervals = pairwise_intervals(rec_pts)

    # global Lorentz boost of every located chart datum
    lam = boost_matrix(0.6, axis=1)
    boosted = [
        LocatedRecord(r.token, r.stage, 0.0, np.zeros(4), r.radius, r.certified)
        for r in records
    ]
    for src, dst in zip(records, boosted):
        chart = np.concatenate(([src.tau], src.x_e[1:]))
        new = lam @ chart
        dst.tau = float(new[0])
        dst.x_e = np.concatenate(([src.x_e[0]], new[1:]))
    classes_b, _ = event_classes(boosted)
    assert sorted(map(sorted, classes_b)) == sorted(map(sorted, classes))
    pts_b = class_positions(boosted, classes_b)
    intervals_b = pairwise_intervals(pts_b)
    assert np.allclose(intervals, intervals_b, atol=0.05, rtol=0.05)

    # gauge relabeling of tokens
    n = len(events)
    perm = {t: (t * 3) % n for t in range(n)}
    relabeled = [
        LocatedRecord(perm[r.token], r.stage, r.tau, r.x_e, r.radius, r.certified)
        for r in records
    ]
    classes_g, _ = event_classes(relabeled)
    expected = sorted(sorted(perm[t] for t in cls) for cls in classes)
    assert sorted(map(sorted, classes_g)) == expected

    # schedule permutation: record processing order is irrelevant
    reordered = list(reversed(records))
    classes_s, _ = event_classes(reordered)
    assert sorted(map(sorted, classes_s)) == sorted(map(sorted, classes))


# ---------------------------------------------------------------------------
# countermodels (prop:event-countermodels)
# ---------------------------------------------------------------------------

def test_population_deficit_drops_event_dimension():
    events, taus = worldline_events(24)
    records = generate_records(events, taus, STAGES, seed=51)
    classes, _ = event_classes(records)
    pts = class_positions(records, classes)
    assert pca_dimension(pts, rel_threshold=5e-2) == 1


def test_label_inflation_doubles_event_count():
    events, taus = generic_events(6, seed=61)
    plain = generate_records(events, taus, STAGES, seed=62)
    classes_plain, _ = event_classes(plain)
    inflated = label_inflated_records(events, taus, STAGES, seed=62)
    classes_inflated, ambiguous = event_classes(inflated)
    assert len(classes_inflated) == 2 * len(classes_plain)
    assert not ambiguous


def test_reconciliation_branching_is_flagged_ambiguous():
    records = branching_records(STAGES)
    germs = germs_from_records(records)
    verdict = pair_verdict(germs[0], germs[1])
    assert verdict == "AMBIGUOUS"
    classes, ambiguous = event_classes(records)
    assert ambiguous == [(0, 1)]
