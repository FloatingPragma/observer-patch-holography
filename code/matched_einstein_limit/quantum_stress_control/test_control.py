from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import measure


def test_optical_commutator_on_entire_occupation_subspace():
    vectors = {"s": np.array([0.1 + 0.2j, 0.5 - 0.7j]),
               "r": np.array([-0.4 + 0.3j, 0.9 + 0.2j])}
    _, levels, bank = measure.maxwell.coordinate_bank(vectors)
    response = 1j * (bank["r"] @ bank["s"] - bank["s"] @ bank["r"])
    expected = -2 * np.vdot(vectors["r"], vectors["s"]).imag
    eye = np.eye(len(levels))
    inside = levels < measure.maxwell.CAP
    np.testing.assert_allclose(response.toarray()[:, inside], expected * eye[:, inside], atol=1e-12)
    # The upper boundary makes the same assertion false if extended to all occupations.
    assert np.linalg.norm((response.toarray() - expected * eye)[:, ~inside]) > 0.1


def test_pinned_packet_stress_and_receiver_draws():
    here = Path(measure.__file__).parent
    packet = json.loads((here / "results.json").read_text())
    assert packet["protocol_sha256"] == sha256((here / "protocol.json").read_bytes()).hexdigest()
    assert packet["parent_code_sha256"] == sha256(Path(measure.maxwell.__file__).read_bytes()).hexdigest()
    n_square = measure.PROTOCOL["enrollments_per_quadrature"]
    n_kick = measure.PROTOCOL["enrollments_per_kick_sign"]
    for run in packet["runs"]:
        assert len(run["weak_optical_channels"]) == 24
        assert run["weak_response_maximum_difference"] < 1e-12
        assert run["maximum_clock_control_difference"] > 0.1
        for stress in run["stress_records"]:
            squares = np.array([[r["sum_of_squared_outcomes_vacuum_photon"] for r in family]
                                for family in stress["field_square_groups_electric_magnetic"]])
            excess = np.diff(squares, axis=-1).squeeze(-1).sum(axis=0) / n_square
            trace = excess.sum() / 2
            np.testing.assert_allclose(trace - excess, stress["record_diagonal_stress"], atol=1e-12)
            np.testing.assert_allclose(sum(stress["normal_ordered_diagonal_stress"]),
                                       run["normal_ordered_energy"], atol=1e-10)
            assert abs(trace - run["normal_ordered_energy"]) < 5 * stress["trace_standard_error"]
        kicks = run["adjacent_xz_finite_kicks_at_half_side_time"]
        gap = []
        for kick in kicks:
            counts = np.array(kick["counts_state_by_kick_sign_by_outcome"])
            np.testing.assert_array_equal(counts.sum(axis=-1), np.full((2, 2), n_kick))
            means = counts @ np.array(kick["receiver_outcome_values"]) / n_kick
            measured = (means[:, 0] - means[:, 1]) / (2 * kick["kick"])
            np.testing.assert_allclose(measured, kick["record_slopes_vacuum_photon"], atol=1e-12)
            gap.append(abs(np.diff(kick["slopes_vacuum_photon"]).item()))
        assert all(12 < a / b < 20 for a, b in zip(gap[:-1], gap[1:], strict=True))
