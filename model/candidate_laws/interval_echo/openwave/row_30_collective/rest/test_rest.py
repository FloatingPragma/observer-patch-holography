"""Finite ring sign and zero-momentum readout checks."""

from pathlib import Path

import numpy as np

from .physics import full_spin_check, propagation
from .read import decode

HERE = Path(__file__).resolve().parent


def test_closing_link_and_physical_sector():
    for n in (5, 6):
        for flip in (-1, 1):
            assert full_spin_check(n, flip, 8.0) < 1e-11
            u, e = propagation(n, 8.0, flip)
            uniform = np.full(n, 1 / np.sqrt(n))
            residual = np.linalg.norm(u @ ((e**2 - 15**2) * (u.T @ uniform)))
            if flip == -((-1) ** n):
                assert residual < 1e-11
            else:
                assert residual > 1


def test_detector_resolves_zero_wave_number_on_both_rings():
    for seed in (312142, 927016):
        with np.load(HERE / f"counts_{seed}.npz", allow_pickle=False) as data:
            a, b = (decode(data, arm) for arm in ("native", "changed"))
        assert 14.996789 < a["rest"] < 15.004576
        assert abs(b["rest"] - 7) < 0.0048
        assert a["spatial"][0]["uniform_weight"] > 0.999
        assert max(m["uniform_weight"] for m in a["spatial"][1:]) < 0.001
        assert max(a["held_sum_error"], a["max_difference_line_error"]) < 0.0048
