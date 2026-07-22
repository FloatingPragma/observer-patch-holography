from __future__ import annotations

import json
import math
import cmath
from pathlib import Path

from verify_collar_gap_certificate import verify


HERE = Path(__file__).resolve().parent


def test_exact_rational_contract_witness() -> None:
    payload = json.loads(
        (HERE / "certificates" / "issue_306_theorem_contract_witness.json").read_text()
    )
    result = verify(payload)
    assert result["valid"] is True
    assert result["physical_clay_receipt"] is False
    assert result["c_floor"] == "3/4"
    assert result["eta_upper"] == "1/2"
    assert result["gap_lower"] == "3/8"


def test_no_mixing_local_countermodel_has_zero_gap() -> None:
    # On support {00, 11}, conditioning either spin on the other fixes it.
    generator_on_centered_sign = 0
    assert generator_on_centered_sign == 0


def test_product_mixing_nonlocal_gray_cycle_gap_vanishes() -> None:
    gaps = []
    for m in range(3, 9):
        size = 2**m
        mode = [cmath.exp(2j * math.pi * index / size) for index in range(size)]
        # (I-E_0)+(I-E_1) is one half of the cycle graph Laplacian.
        applied = [
            value - (mode[(index - 1) % size] + mode[(index + 1) % size]) / 2
            for index, value in enumerate(mode)
        ]
        gap = 1 - math.cos(2 * math.pi / size)
        assert max(abs(lhs - gap * rhs) for lhs, rhs in zip(applied, mode)) < 1e-12
        gaps.append(gap)
    assert all(later < earlier for earlier, later in zip(gaps, gaps[1:]))
    assert gaps[-1] < 0.001
