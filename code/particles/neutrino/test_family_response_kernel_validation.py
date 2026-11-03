"""Validation controls for the supplied symmetric neutrino kernel.

The family-response producer used to replace every kernel diagonal entry below
``1e-15`` by ``1e-15`` and to clip the normalized matrix before taking
``-log``.  That silently accepted non-PSD kernels and corrupted valid kernels
with small entries.  These controls use the producer's command line so they
exercise the same path as the committed artifacts.
"""

from __future__ import annotations

import json
from decimal import Decimal, localcontext
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FAMILY_PRODUCER = ROOT / "particles" / "neutrino" / "derive_family_response_tensor.py"
HOLONOMY_PRODUCER = ROOT / "particles" / "neutrino" / "derive_majorana_holonomy_lift.py"
SECTOR = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
COMMITTED_FAMILY = ROOT / "particles" / "runs" / "neutrino" / "family_response_tensor.json"


def _sector_with_kernel(tmp_path: Path, kernel: list[list[float]]) -> Path:
    payload = json.loads(SECTOR.read_text(encoding="utf-8"))
    nu = payload["sector_response_object"]["nu"]
    for key in ("K_core", "K_core_directed", "K_core_majorana_sym"):
        nu[key] = kernel
    path = tmp_path / "sector.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _run_family(tmp_path: Path, kernel: list[list[float]]):
    sector = _sector_with_kernel(tmp_path, kernel)
    output = tmp_path / "family_response.json"
    result = subprocess.run(
        [sys.executable, str(FAMILY_PRODUCER), "--sector", str(sector), "--out", str(output)],
        capture_output=True, text=True,
    )
    return result, output


def test_nonpositive_diagonal_kernel_is_rejected(tmp_path: Path) -> None:
    kernel = [[1.0, 0.1, 0.1], [0.1, -1.0, 0.1], [0.1, 0.1, 1.0]]
    result, output = _run_family(tmp_path, kernel)
    assert result.returncode != 0
    assert "strictly positive diagonal" in result.stderr
    assert not output.exists()


def test_positive_diagonal_non_psd_kernel_is_rejected(tmp_path: Path) -> None:
    # Symmetric with positive diagonal, eigenvalue 1 - 2 = -1.
    kernel = [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    result, output = _run_family(tmp_path, kernel)
    assert result.returncode != 0
    assert "positive semidefinite" in result.stderr
    assert not output.exists()


def test_nonpositive_normalized_entry_is_rejected(tmp_path: Path) -> None:
    # Positive-definite correlation kernel whose 12 entry is negative: the
    # -log defect weight is undefined there, and no floor may replace it.
    kernel = [[1.0, -0.5, 0.1], [-0.5, 1.0, 0.1], [0.1, 0.1, 1.0]]
    result, output = _run_family(tmp_path, kernel)
    assert result.returncode != 0
    assert "strictly positive entries" in result.stderr
    assert not output.exists()


def test_holonomy_producer_rejects_nonpositive_diagonal(tmp_path: Path) -> None:
    sector = _sector_with_kernel(tmp_path, [[1.0, 0.1, 0.1], [0.1, -1.0, 0.1], [0.1, 0.1, 1.0]])
    output = tmp_path / "holonomy.json"
    result = subprocess.run(
        [sys.executable, str(HOLONOMY_PRODUCER), "--input", str(sector), "--output", str(output)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "strictly positive diagonal" in result.stderr
    assert not output.exists()


def test_valid_small_kernel_is_normalized_exactly_without_a_floor(tmp_path: Path) -> None:
    """Independent rational control for a kernel below the old 1e-15 floor."""
    correlations = [[1, 0.25, 0.125], [0.25, 1, 0.5], [0.125, 0.5, 1]]
    scale = [1e-20, 1e-19, 1e-18]
    kernel = [[scale[i] * scale[j] * correlations[i][j] for j in range(3)] for i in range(3)]
    result, output = _run_family(tmp_path, kernel)
    assert result.returncode == 0, result.stderr
    payload = json.loads(output.read_text(encoding="utf-8"))

    np.testing.assert_allclose(payload["majorana_normalization_diag"], scale, rtol=1e-15)
    np.testing.assert_allclose(payload["E_nu"], correlations, rtol=1e-12, atol=0.0)

    with localcontext() as context:
        context.prec = 60
        expected_weights = [
            [-Decimal(correlations[i][j]).ln() for j in range(3)] for i in range(3)
        ]
    np.testing.assert_allclose(
        payload["S_nu"], [[float(value) for value in row] for row in expected_weights],
        rtol=1e-12, atol=1e-14,
    )
    # The old floor returned 3.16e-8 here and defect weights near 57 instead.
    assert max(payload["majorana_normalization_diag"]) < 1e-17


def test_committed_artifact_is_unchanged(tmp_path: Path) -> None:
    output = tmp_path / "family_response.json"
    result = subprocess.run(
        [sys.executable, str(FAMILY_PRODUCER), "--out", str(output)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    fresh = json.loads(output.read_text(encoding="utf-8"))
    committed = json.loads(COMMITTED_FAMILY.read_text(encoding="utf-8"))
    fresh.pop("generated_utc")
    committed.pop("generated_utc")
    assert committed == fresh
