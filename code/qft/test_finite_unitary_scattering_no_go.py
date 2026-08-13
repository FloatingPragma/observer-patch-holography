from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

from finite_unitary_scattering_no_go import (
    audit_step,
    default_replays,
    identity,
    relative_comparison,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
LEAN_ROOT = REPO_ROOT / "Lean"
LEAN_SOURCE = LEAN_ROOT / "QFT" / "FiniteUnitaryScatteringNoGo.lean"
BOUNDARY_SURFACE = (
    REPO_ROOT / "paper" / "tex_fragments" / "QFT_STRUCTURAL_INHERITANCE_STATUS.tex"
)


def test_exact_controls_attack_both_claim_and_boundary() -> None:
    identity_replay, swap, quarter_turn = default_replays(horizon=17)

    assert not identity_replay.nonidentity
    assert identity_replay.adjacent_displacements_squared == (0,) * 17

    for replay in (swap, quarter_turn):
        assert replay.nonidentity
        assert replay.one_step_displacement_squared > 0
        assert replay.adjacent_displacements_squared == (
            replay.one_step_displacement_squared,
        ) * 17
        assert replay.factorization_holds
        assert replay.identical_relative_comparison_is_constant


def test_relative_comparison_control_can_converge_for_nontrivial_steps() -> None:
    swap = ((0, 1), (1, 0))
    for exponent in range(25):
        assert relative_comparison(swap, swap, exponent) == identity(2)


def test_nonorthogonal_input_fails_closed() -> None:
    with pytest.raises(ValueError, match="exact orthogonal"):
        audit_step("shear", ((1, 1), (0, 1)))


def test_noninteger_input_fails_closed_instead_of_truncating() -> None:
    with pytest.raises(ValueError, match="exact integers"):
        audit_step("fractional", ((1, 0.5), (0, 1)))


def test_lean_source_has_no_proof_escape_hatches() -> None:
    source = LEAN_SOURCE.read_text(encoding="utf-8")
    forbidden = re.compile(r"(?m)^\s*(?:axiom|admit)\b|\bsorry\b")
    assert not forbidden.search(source)
    assert "finite_unitary_ambient_powers_have_no_limit" in source
    assert "identical_relative_evolution_tendsto" in source


def test_paper_boundary_forbids_scattering_promotion() -> None:
    boundary = " ".join(BOUNDARY_SURFACE.read_text(encoding="utf-8").split())
    assert "relative evolution can converge" in boundary
    assert "Routes not excluded here include comparison dynamics" in boundary
    assert "proves no S-matrix" in boundary
    assert r"it does not close issue~\#743" in boundary


def test_lean_kernel_checks_the_result() -> None:
    lake = shutil.which("lake")
    if lake is None:
        pytest.skip("Lean lake executable is unavailable")
    completed = subprocess.run(
        [lake, "env", "lean", "QFT/FiniteUnitaryScatteringNoGo.lean"],
        cwd=LEAN_ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "sorryAx" not in completed.stdout + completed.stderr
