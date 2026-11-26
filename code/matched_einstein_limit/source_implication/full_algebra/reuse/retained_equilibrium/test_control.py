"""Independent capped records, native gates, return memory and held response."""

from itertools import product
import json

import numpy as np
import pytest

from . import certificate, control, gates, held, measure, probe


@pytest.fixture(scope="module")
def design():
    return control.compare()


def test_all_24_terminal_bits_against_floating_propagation():
    bank = certificate.channels()
    x, *_ = probe.density_center()
    for iq, p in product(range(2), repeat=2):
        raw = probe.wave_bank((-1, 1)[iq], bool(p))
        for j, route in enumerate(probe.reuse.ROUTES):
            y = np.zeros((2, 4, 4), complex)
            for (kind, a, b), ks in probe.histories(raw, route):
                y[int(kind == "success" and a + b < 12)] += probe.conditional(ks, x[p])
            result = control.state_outputs(bank[iq, p, j], x[p])
            np.testing.assert_allclose(y, result, atol=5e-10, rtol=0)


def test_majorization_and_mixing(design):
    x, schemes, rows = design
    for p, row in enumerate(rows):
        assert min(row["with_bit"]) > 0.00639 and row["without_bit"][1] < -0.00465
        assert row["bit_only"][1] < -0.00748 and row["maximum_uniform_twirl"] > 0.1747
        slots = schemes[p].sum(axis=0)
        np.testing.assert_allclose(
            slots @ x[p].ravel(), x[p].ravel(), atol=3e-12, rtol=0
        )
        check = slots @ np.eye(4).ravel() * x[p].trace() / 4
        for _ in range(25):
            check = slots @ check
        assert np.linalg.norm(check - x[p].ravel()) < 1e-5
        assert row["maximum_gate_error"] < 2e-12
    ans = measure.successive(schemes, x)
    assert (
        ans["contrast_band_99"][0]
        < ans["contrast"]
        < ans["contrast_band_99"][1]
        < -0.00013
    )


def test_shuttle_and_native_rotations():
    assert [len(p) - 1 for p in gates.shuttle()] == [18, 16]
    cz = gates.native_cz()
    np.testing.assert_allclose(cz, np.diag([1, 1, 1, -1]), atol=1e-13)
    rng = np.random.default_rng(2026)
    for _ in range(5):
        u = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[0]
        np.testing.assert_allclose(gates.compile_gate(u, cz), u, atol=1e-12, rtol=0)


def test_held_balanced_stationary_sources(design):
    x, ch, _ = design
    observed = held.read(ch, x)
    archive = json.loads((certificate.HERE / "results.json").read_text())["held"]
    assert observed["fixed_state_trace_error"] < 3e-10
    assert observed["phase_axes"] == [2, 0]
    assert abs(observed["volume"]) < 1e-9
    assert 0.4062 < observed["entropy"] < 0.4063
    assert abs(observed["linear_residual"]) < 1e-6
    for actual, saved in zip(observed["finite"], archive["finite"], strict=True):
        np.testing.assert_allclose(
            [actual["area"], actual["entropy"]],
            [saved["area"], saved["entropy"]],
            atol=1e-5,
            rtol=0,
        )
    assert (
        abs(observed["finite"][1]["residual"])
        < abs(observed["finite"][0]["residual"]) / 4
    )
