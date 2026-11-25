import json

import numpy as np

from . import live_source as live
from .source import measure as s


def test_same_marks_retained_with_inclusive_trace_and_coupling_control():
    row = next(
        r
        for r in json.loads((s.HERE / "results.json").read_text())["runs"]
        if (r["level"], r["time"], r["forming_charge"], r["pump_charge"])
        == (1, 4, -1, -1)
    )
    x, _, z = s.formation.PAULI
    for coupled in (True, False):
        kraus = live.head_channel(coupled)
        assert kraus.shape == (960, 8, 8)
        np.testing.assert_allclose(
            np.einsum("pai,paj->ij", kraus.conj(), kraus), np.eye(8), rtol=0, atol=2e-12
        )
        states = live.postpump(coupled)
        for name, rho in states.items():
            assert abs(np.trace(rho) - (name != "hidden")) < 2e-12
            np.testing.assert_allclose(rho, rho.conj().T, atol=1e-14, rtol=0)
            tensor = rho.reshape(2, 4, 2, 4)
            strain = np.trace(tensor, axis1=1, axis2=3)
            heads = np.trace(tensor, axis1=0, axis2=2)
            if name != "hidden":
                assert np.min(np.linalg.eigvalsh(rho)) >= -2e-12
            if coupled and name != "hidden":
                archived = row["arms"][name]
                expected = np.array(archived["density_real"]) + 1j * np.array(
                    archived["density_imag"]
                )
                np.testing.assert_allclose(strain, expected, atol=2e-12, rtol=0)
            elif not coupled:
                np.testing.assert_allclose(
                    rho, np.kron(np.ones((2, 2)) / 2, heads), atol=2e-12, rtol=0
                )
            else:
                np.testing.assert_allclose(strain, 0, atol=2e-12, rtol=0)
                assert -0.11579 < np.trace(heads @ np.kron(x, z)).real < -0.11577
        # Replacing input heads by independent copies preserves both marginals;
        # after the pump the output joint state records strain–head correlations.
        joint = states["reference"].reshape(2, 4, 2, 4)
        separated = np.kron(
            np.trace(joint, axis1=1, axis2=3), np.trace(joint, axis1=0, axis2=2)
        )
        if coupled:
            assert np.linalg.norm(states["reference"] - separated) > 0.01


def test_normalized_inverse_witness_and_certified_live_densities():
    data = live.certify()
    assert data["kernel_basis"] == ["ZZ", "XZ", "ZY", "ZX"]
    assert data["certified_two_sided_step"] == 1e-8
    assert data["normalized_direction_l2_error_upper"] < 1.35e-9
    assert data["refined_unit_center"]["l2_error"] < 1e-13
    assert (
        0.03401
        < data["hidden_modular_nats"][0]
        < data["hidden_modular_nats"][1]
        < 0.03403
    )
    for control in data["channels"].values():
        for arm, err in control["nominal_state_trace_norm_errors"].items():
            assert err < (1.35e-9 if arm == "hidden" else 3e-12)
        assert control["kraus_channel_trace_norm_error_per_unit_input"] < 1e-9
    calibration = live.source_calibration()
    np.testing.assert_allclose(
        calibration["ordinary_modular"], np.log(2) / 3, rtol=0, atol=1e-15
    )
    assert (
        data["decorrelation_modular_nats"][0]
        < calibration["decorrelated_modular"]
        < data["decorrelation_modular_nats"][1]
    )
