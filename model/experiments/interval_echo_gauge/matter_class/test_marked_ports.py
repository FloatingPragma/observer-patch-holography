"""Prospective marked-current forecasts from the frozen B14 port action."""

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import runpy

import numpy as np


def test_marked_port_ramsey_forecasts():
    root = Path(__file__).resolve().parents[4]
    b14 = root / "code/b14_jacobi"
    selector = runpy.run_path(str(b14 / "port_dual_metric_selector.py"))
    cert = json.loads((b14 / "port_dual_metric_selector.certificate.json").read_text())
    receipt = selector["RECEIPT"]
    assert (
        sha256(receipt.read_bytes()).hexdigest() == cert["pinned_inputs"][receipt.name]
    )
    projectors = selector["spectral_projectors"](selector["group_rows"]())
    Q5 = selector["Q5"]

    # Reference basis is (positive port-0 spin, orthogonal spin, singlet).
    spin_z = np.diag([1, -1])
    spin_x = np.array([[0, 1], [1, 0]])
    bright = np.array([1, 0, 1]) / np.sqrt(2)
    dark = np.array([-1, 0, 1]) / np.sqrt(2)
    analyzer = np.array([bright, dark, [0, 1, 0]])
    for sector, sign in (("three_plus", 1), ("three_minus", -1)):
        p = projectors[sector]
        assert p[0][0] == p[1][1] == Q5(Fraction(1, 4))
        assert p[0][1] == Q5(0, Fraction(sign, 20))
        c = p[0][1] / p[0][0]
        assert c == Q5(0, Fraction(sign, 5))
        # The cross product phase is fixed by chronological matrix multiplication.
        cosine = sign / np.sqrt(5)
        axis1 = cosine * spin_z + np.sqrt(1 - cosine**2) * spin_x
        for polarity in (1, -1):
            full = np.eye(3, dtype=complex)
            full[:2, :2] = (-1j * polarity * axis1) @ (-1j * spin_z)
            probs = abs(analyzer @ full @ bright) ** 2
            expected = [
                Q5(Fraction(3, 10)) - polarity * c / 2,
                Q5(Fraction(3, 10)) + polarity * c / 2,
                Q5(Fraction(2, 5)),
            ]
            assert sum(expected) == Q5(1)
            assert expected[0] - expected[1] == -polarity * c
            numeric = [float(t.a) + np.sqrt(5) * float(t.b) for t in expected]
            np.testing.assert_allclose(probs, numeric, atol=1e-14, rtol=0)
