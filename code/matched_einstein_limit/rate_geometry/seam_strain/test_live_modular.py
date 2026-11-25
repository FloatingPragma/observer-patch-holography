import json

import numpy as np
from scipy.linalg import expm

from . import live_modular, live_source, live_survey


def test_retained_modular_operator_and_signed_response_enclosures():
    actual = live_modular.certify()
    archive = json.loads(
        (live_modular.source.HERE.parent / "live_modular.json").read_text()
    )
    ref, deltas, _ = live_survey.retained(True)
    for data in (archive, actual):
        assert data["reference_minimum_eigenvalue_lower"] >= 0.0014
        assert data["reference_log_op_error"] < 5e-11
        k = np.array(data["reference_log_real"]) + 1j * np.array(
            data["reference_log_imag"]
        )
        np.testing.assert_allclose(k, k.conj().T, atol=2e-14, rtol=0)
        np.testing.assert_allclose(expm(-k), ref, atol=2e-12, rtol=0)
        expected = {
            "ordinary": 0.02037788,
            "decorrelated": 0.04089855,
            "hidden": -0.00112805,
            "left": 1.09658e-6,
            "right": -1.18279e-5,
        }
        for name, center in expected.items():
            low, high = data["modular_slopes"][name]
            assert (
                high >= low
                and high - low < 2e-10
                and abs((high + low) / 2 - center) < 1e-8
            )
            if data is actual:
                # Current nominal hidden uses the simple XZ+ZX center. Apply
                # the certified trace error before comparing this finite read.
                allowance = 7 * data["nominal_tangent_trace_norm_errors"][name]
                direction = float(np.trace(k @ deltas[name]).real)
                assert low - allowance - 1e-12 <= direction <= high + allowance + 1e-12
            a, b = archive["modular_slopes"][name]
            assert max(a, low) <= min(b, high)
        for name in ("left", "right"):
            assert data["rotation_trace_norm_errors"][name] < 1e-11
            a, b = data["rotation_output_trace_norms"][name]
            assert 0 < a <= b and b - a < 1e-13
            n = np.linalg.svd(deltas[name], compute_uv=False).sum()
            assert a - 1e-11 <= n <= b + 1e-11


def test_local_head_rotations_preserve_marginals_and_incoming_modular_slope():
    form = live_modular.source.formation
    board = form.federation(1)
    rho, _ = live_modular.source.form(board, form.region(board), -1)
    k = -live_modular.source.pairing.log_density(rho)
    w = live_source.head_channel()
    for h in live_survey.rotations(rho).values():
        d = -1j * (h @ rho - rho @ h)
        assert max(np.linalg.norm(x) for x in form.marginals(d)) < 5e-15
        assert abs(np.trace(k @ d)) < 1e-13
        # Check the finite orbit, not only the incoming commutator formula.
        t = 0.25
        u = expm(-1j * t * h)
        rotated = u @ rho @ u.conj().T
        for p, q in zip(form.marginals(rotated), form.marginals(rho)):
            np.testing.assert_allclose(p, q, atol=2e-14, rtol=0)
        assert np.linalg.norm(live_source.apply_channel(w, d, (2 / 3, 1 / 3))) > 1e-5
