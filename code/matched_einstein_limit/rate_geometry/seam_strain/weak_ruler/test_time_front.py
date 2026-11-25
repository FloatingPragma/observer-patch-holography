import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import flight, late_read, observability, time_front
from .. import live_source
from ..source import measure as source


def test_front_ranks_replay_from_the_hamiltonian():
    data = time_front.run()
    saved = json.loads((flight.HERE / "time_front.json").read_text())
    for name in ("with_front_derivatives", "with_delayed_fronts"):
        assert data[name]["rank"] == saved[name]["rank"] == 9
        assert (
            data[name]["inverse_defect_upper"]
            == saved[name]["inverse_defect_upper"]
            < 0.01
        )
    for record in (data, saved):
        for lo, hi in record["matched_age_four_lotteries"]["delayed_front_contrasts"]:
            assert lo < hi and lo * hi > 0
        witness = record["matched_preparations"]
        assert witness["rank"] == 5
        assert witness["fixed_readings"] == [
            "z",
            "dz",
            "ddz",
            "outward_number_radius_1",
            "outward_number_radius_2",
        ]
        assert all(lo > 0 for lo, _ in witness["signed_responses"].values())
    _, base, _, _ = time_front.measurements()
    m, e = time_front.stack(base)
    independent = late_read.front_response_matrix()
    _, pump_error = late_read.exact_pump()
    assert all(
        abs(a - b) <= error + pump_error
        for a, b, error in zip(m, independent, e, strict=True)
    )


def test_native_pump_clock_and_finite_differences_on_a_held_preparation():
    # Propagate the source at three pump ages on an independent floating-point
    # Krylov solver. The optical operators and launch preparation are shared.
    board = source.formation.federation(1)
    formed, product = source.form(board, source.formation.region(board), -1)
    head_change = product - formed
    inputs = np.kron(np.diag([2 / 3, 1 / 3]), head_change)
    coordinates = np.array(
        [
            np.trace(4 * np.asarray(axis, complex) @ head_change).real
            for axis in observability.AXES
        ]
    )
    h = -1j * flight.hardware(-1)[0] / 40
    kraus = live_source.head_channel(True)
    wave = kraus.reshape(-1, 2, 4, 8).transpose(1, 0, 2, 3).reshape(-1, 8)

    def observed(dt, ops):
        y = expm_multiply(dt * h, wave, traceA=dt * h.diagonal().sum()) if dt else wave
        blocks = y.reshape(2, -1, 4, 8).transpose(1, 0, 2, 3).reshape(-1, 8, 8)
        rho = np.einsum("pai,ij,pbj->ab", blocks, inputs, blocks.conj(), optimize=True)
        return np.array([np.trace(op @ rho).real for op in ops])

    local, _, rates, delayed = time_front.measurements()
    operators = time_front.operators()
    selected = [
        np.array(operators[i][0], float) / float(operators[i][1])
        for i in time_front.REPEATED
    ]
    probes = [np.diag([1] * 4 + [-1] * 4), *selected]
    dt = 0.0001
    left, mid, right = [observed(t, probes) for t in (-dt, 0, dt)]
    first = (right - left) / (2 * dt)
    acceleration = (left[0] + right[0] - 2 * mid[0]) / dt**2

    def prediction(rows):
        return np.array(time_front.stack(rows)[0], float) @ coordinates

    np.testing.assert_allclose(
        prediction(local), [mid[0], first[0], acceleration], rtol=5e-5, atol=8e-8
    )
    np.testing.assert_allclose(prediction(rates), first[1:], rtol=2e-5, atol=2e-10)
    np.testing.assert_allclose(
        prediction(delayed), observed(0.125, selected), rtol=1e-10, atol=3e-13
    )
