#!/usr/bin/env python3
"""Regression tests for GitHub issue #530's triangle-free cycle counterexample."""

from __future__ import annotations

import re
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

from transport_from_gluing import (  # noqa: E402
    Matrix,
    identity_matrix,
    multiply,
    path_transport,
    strict_transportability_report,
)


IDENTITY_2 = identity_matrix(2)
NONTRIVIAL_HOLONOMY: Matrix = (
    (-1.0 + 0.0j, 0.0 + 0.0j),
    (0.0 + 0.0j, 1.0 + 0.0j),
)


def triangle_free_cycle() -> tuple[dict[tuple[int, int], Matrix], tuple[tuple[int, int, int], ...]]:
    transports = {
        (0, 1): IDENTITY_2,
        (1, 0): IDENTITY_2,
        (1, 2): IDENTITY_2,
        (2, 1): IDENTITY_2,
        (2, 3): IDENTITY_2,
        (3, 2): IDENTITY_2,
        (3, 0): NONTRIVIAL_HOLONOMY,
        (0, 3): NONTRIVIAL_HOLONOMY,
    }
    return transports, ()


def test_triangle_free_cycle_rejects_nontrivial_holonomy_after_central_strictification() -> None:
    transports, two_simplices = triangle_free_cycle()
    central_degree_two_class_trivial = len(two_simplices) == 0
    cycle_holonomy = path_transport((0, 1, 2, 3, 0), transports)

    report = strict_transportability_report(
        associator_defect_strictifiable=central_degree_two_class_trivial,
        fundamental_loop_holonomies=(cycle_holonomy,),
        expected_loop_generator_count=1,
    )

    assert report.associator_defect_strictifiable is True
    assert report.represented_loop_holonomy_trivial is False
    assert report.strict_path_independent is False
    assert cycle_holonomy[0][0] != cycle_holonomy[1][1]  # no central rephasing removes it
    assert path_transport((0, 1, 2), transports) != path_transport((0, 3, 2), transports)


def test_triangle_free_cycle_rejects_nontrivial_holonomy_after_noncentral_strictification() -> None:
    transports, two_simplices = triangle_free_cycle()
    higher_associator_strictifiable = len(two_simplices) == 0
    strict_g_cocycle_holonomy = path_transport((0, 1, 2, 3, 0), transports)
    band_holonomy = (
        strict_g_cocycle_holonomy[0][0] * strict_g_cocycle_holonomy[1][1]
        - strict_g_cocycle_holonomy[0][1] * strict_g_cocycle_holonomy[1][0]
    )

    report = strict_transportability_report(
        associator_defect_strictifiable=higher_associator_strictifiable,
        fundamental_loop_holonomies=(strict_g_cocycle_holonomy,),
        expected_loop_generator_count=1,
    )

    assert report.associator_defect_strictifiable is True
    assert report.represented_loop_holonomy_trivial is False
    assert report.strict_path_independent is False
    # For SU(2) -> U(2), determinant detects the U(1) band and is unchanged by SU(2) edge moves.
    assert band_holonomy == -1.0 + 0.0j


def test_path_transport_uses_manuscript_order_for_noncommuting_edges() -> None:
    pauli_x: Matrix = ((0.0j, 1.0 + 0.0j), (1.0 + 0.0j, 0.0j))
    pauli_z: Matrix = ((1.0 + 0.0j, 0.0j), (0.0j, -1.0 + 0.0j))
    transports = {(0, 1): pauli_x, (1, 2): pauli_z}

    composite = path_transport((0, 1, 2), transports)

    assert composite == multiply(pauli_z, pauli_x)
    assert composite != multiply(pauli_x, pauli_z)


def _collapsed(path: Path) -> str:
    return re.sub(r"\s+", " ", path.read_text(encoding="utf-8"))


def test_paper_sources_retain_the_combined_criterion_and_cycle_check() -> None:
    gauge_wrapper = _collapsed(
        ROOT
        / "paper"
        / "deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex"
    )
    shared_spine = _collapsed(ROOT / "paper" / "tex_fragments" / "PAPER.tex")
    gauge = f"{gauge_wrapper} {shared_spine}"

    assert "\\input{tex_fragments/PAPER.tex}" in gauge_wrapper
    assert "\\def\\OPHSkipEinsteinChain{1}" in gauge_wrapper
    assert "Triangle-free cycle check" in gauge
    assert "associator-strictification and trivial-holonomy criterion" in gauge
    assert "full orbit \\(q_\\Sigma\\)" in gauge
    assert "\\(o^{(2)}_\\Sigma=0\\)" in gauge
    assert "\\SU(2)\\hookrightarrow\\U(2)" in gauge
    assert "fixed abelian unitary coefficient group \\(Z_\\Sigma\\)" in gauge
    assert "central defect strictification and residual holonomy" in gauge

    assert "strict path-independent transport exists if and only if \\([z]_\\Sigma=0\\)" not in gauge
    assert "transport of collar charges exists if and only if \\([z]_\\Sigma=0\\)" not in gauge
    assert "strict path-independent transport exists iff the central loop-coherence class" not in gauge
    assert "Strict ordinary transport exists iff \\(q_\\Sigma=[(g,h)]\\) vanishes." not in gauge
    assert "\\(q_\\Sigma=0\\) strictifies the higher defect" not in gauge
    assert "A loop-coherent global gluing exists iff {[}z{]} = 0" not in gauge
