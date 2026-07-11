#!/usr/bin/env python3
"""Tests for the quotient-intrinsic geometry producer receipts (GitHub #523).

Each test verifies one clause of the compact paper's
`subsec:quotient-intrinsic-geometry-producer`.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from quotient_cap_readout import (  # noqa: E402
    IncidenceComplex,
    RepairSystem,
    boundary_4_simplex_2_skeleton,
    classify_surface,
    complexes_equal,
    complexes_isomorphic_under,
    cross_ratio,
    csaszar_torus,
    euler_characteristic,
    gauge_relabel,
    icosahedron,
    incidence_complex,
    is_closed_surface,
    is_connected,
    kms_receipt,
    minkowski,
    mobius_normalize,
    normal_form_records,
    orient,
    produced_cap_normal,
    readout_from_system,
    reconstruct_from_cross_ratios,
    refinement_is_simplicial,
    refinement_subdivide,
    spherical_incidence_receipt,
    stereographic,
)


def build_systems() -> dict[str, RepairSystem]:
    ico, _ = icosahedron()
    return {
        "S2": RepairSystem(ico),
        "T2": RepairSystem(csaszar_torus()),
        "S3": RepairSystem(boundary_4_simplex_2_skeleton()),
        "WEDGE": RepairSystem(
            [tuple(t) for t in __import__("quotient_cap_readout").wedge_of_two_spheres()]
        ),
    }


# ---------------------------------------------------------------------------
# D1 clauses: termination, confluence, schedule independence
# ---------------------------------------------------------------------------

def test_repair_terminates_and_is_schedule_independent():
    rng = random.Random(0)
    for name, system in build_systems().items():
        reference = system.repair()
        assert not RepairSystem.conflicted_records(
            type("S", (), {"records": system.records, "state": reference})()
        ), name
        for _ in range(10):
            schedule = list(range(len(system.records)))
            rng.shuffle(schedule)
            assert system.repair(schedule) == reference, name


def test_rewrite_systems_isomorphic_across_countermodels():
    sigs = {name: s.rewrite_signature() for name, s in build_systems().items()}
    assert len(set(sigs.values())) == 1, sigs


# ---------------------------------------------------------------------------
# incidence complex invariance (Lemma lem:incidence-invariance)
# ---------------------------------------------------------------------------

def test_incidence_complex_schedule_invariance():
    system = build_systems()["S2"]
    ref = readout_from_system(system)
    rng = random.Random(1)
    for _ in range(5):
        schedule = list(range(len(system.records)))
        rng.shuffle(schedule)
        assert complexes_equal(ref, readout_from_system(system, schedule))


def test_incidence_complex_gauge_invariance():
    ico, _ = icosahedron()
    system = RepairSystem(ico)
    nf = system.repair()
    k_ref = incidence_complex(normal_form_records(system, nf))
    rng = random.Random(2)
    labels = list(k_ref.vertices)
    shuffled = labels[:]
    rng.shuffle(shuffled)
    perm = dict(zip(labels, shuffled))
    relabeled = gauge_relabel(normal_form_records(system, nf), perm)
    k_gauge = incidence_complex(relabeled)
    assert complexes_isomorphic_under(k_ref, k_gauge, perm)


def test_refinement_naturality():
    ico, _ = icosahedron()
    system = RepairSystem(ico)
    nf = system.repair()
    coarse_records = normal_form_records(system, nf)
    fine_records, projection = refinement_subdivide(coarse_records)
    coarse = incidence_complex(coarse_records)
    fine = incidence_complex(fine_records)
    assert refinement_is_simplicial(fine, coarse, projection)
    # the refined complex still carries the spherical receipt
    assert spherical_incidence_receipt(fine)


# ---------------------------------------------------------------------------
# receipts and topology production (thm:topology-production)
# ---------------------------------------------------------------------------

def test_topology_production_and_underdetermination():
    verdicts = {}
    for name, system in build_systems().items():
        k = readout_from_system(system)
        verdicts[name] = (spherical_incidence_receipt(k), classify_surface(k))
    assert verdicts["S2"] == (True, "S2")
    assert verdicts["T2"] == (False, "T2")
    assert verdicts["S3"][0] is False
    assert verdicts["S3"][1] == "NOT_A_CLOSED_SURFACE"
    assert verdicts["WEDGE"][0] is False
    assert verdicts["WEDGE"][1] == "NOT_A_CLOSED_SURFACE"


def test_euler_characteristics_and_surface_clauses():
    ico, _ = icosahedron()
    k_s2 = incidence_complex(ico)
    assert euler_characteristic(k_s2) == 2
    assert is_closed_surface(k_s2) and is_connected(k_s2)
    assert orient(k_s2) is not None

    k_t2 = incidence_complex(csaszar_torus())
    assert euler_characteristic(k_t2) == 0
    assert is_closed_surface(k_t2) and is_connected(k_t2)
    assert orient(k_t2) is not None

    k_s3 = incidence_complex(boundary_4_simplex_2_skeleton())
    assert not is_closed_surface(k_s3)  # every edge lies in three triangles


def test_wrong_beta_kms_receipt_separates():
    two_pi = 2.0 * np.pi
    assert kms_receipt(two_pi)
    assert not kms_receipt(two_pi * 1.1)
    assert not kms_receipt(1.0)


# ---------------------------------------------------------------------------
# conformal/cap production (thm:conformal-cap-production)
# ---------------------------------------------------------------------------

def test_cross_ratio_reconstruction_of_celestial_embedding():
    _, coords = icosahedron()
    # rotate slightly so no vertex sits at the exact north pole
    theta = 0.3
    rot = np.array(
        [[1, 0, 0],
         [0, np.cos(theta), -np.sin(theta)],
         [0, np.sin(theta), np.cos(theta)]]
    )
    pts = coords @ rot.T
    gauge = (0, 1, 2)
    reconstructed = reconstruct_from_cross_ratios(pts, gauge)
    zs = [stereographic(p) for p in pts]
    g1, g2, g3 = (zs[i] for i in gauge)
    for i, z in enumerate(zs):
        if i == gauge[2]:
            assert np.isinf(reconstructed[i].real)
        else:
            assert abs(reconstructed[i] - mobius_normalize(z, g1, g2, g3)) < 1e-9


def test_cross_ratio_gauge_covariance():
    # cross-ratios are Mobius invariants: applying a Mobius map to all points
    # leaves every receipt cross-ratio unchanged.
    rng = np.random.default_rng(4)
    zs = rng.normal(size=4) + 1j * rng.normal(size=4)

    def moebius(z: complex) -> complex:
        return (2.0 * z + 1.0j) / (1.0j * z + 1.5)

    cr1 = cross_ratio(*zs)
    cr2 = cross_ratio(*[moebius(z) for z in zs])
    assert abs(cr1 - cr2) < 1e-9


def test_produced_cap_normal_matches_paper_formula():
    rng = np.random.default_rng(6)
    c = rng.normal(size=3)
    c /= np.linalg.norm(c)
    alpha = 0.7
    # boundary circle points: c . x = cos(alpha)
    basis = np.linalg.svd(c.reshape(1, 3))[2][1:]
    angles = np.linspace(0.0, 2.0 * np.pi, 9, endpoint=False)
    pts = np.stack(
        [np.cos(alpha) * c
         + np.sin(alpha) * (np.cos(t) * basis[0] + np.sin(t) * basis[1])
         for t in angles]
    )
    n_c = produced_cap_normal(pts)
    expected = np.concatenate(
        ([np.cos(alpha) / np.sin(alpha)], c / np.sin(alpha))
    )
    assert np.allclose(n_c, expected, atol=1e-8)
    # unit spacelike, and sign classification separates interior/exterior
    assert abs(minkowski(n_c, n_c) - 1.0) < 1e-8
    interior = np.concatenate(([1.0], c))
    exterior = np.concatenate(([1.0], -c))
    assert minkowski(n_c, interior) > 0.0
    assert minkowski(n_c, exterior) < 0.0


def test_cap_normal_residual_shrinks_under_refinement():
    rng = np.random.default_rng(7)
    c = np.array([0.0, 0.0, 1.0])
    alpha = 0.7
    basis = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    expected = np.concatenate(([np.cos(alpha) / np.sin(alpha)], c / np.sin(alpha)))
    residuals = []
    for stage, (npts, noise) in enumerate([(6, 1e-2), (12, 1e-3), (24, 1e-4)]):
        angles = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)
        pts = np.stack(
            [np.cos(alpha) * c
             + np.sin(alpha) * (np.cos(t) * basis[0] + np.sin(t) * basis[1])
             for t in angles]
        )
        pts += noise * rng.normal(size=pts.shape)
        pts /= np.linalg.norm(pts, axis=1, keepdims=True)
        residuals.append(np.linalg.norm(produced_cap_normal(pts) - expected))
    assert residuals[0] > residuals[1] > residuals[2]
    assert residuals[2] < 1e-3
