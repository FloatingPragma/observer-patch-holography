import json

import numpy as np
import pytest

from . import certify as c


def test_two_price_charts_and_hidden_entropy_response():
    result = c.certify()
    assert result == json.loads((c.HERE / "results.json").read_text())
    hidden = result["blind_signed_curve"]
    assert hidden["entropy_slope_per_unit"][0] > 6_092_200
    assert hidden["parity_times_yz_over_two_slope"][0] > 1_106_890_000
    assert result["shared_shape_entropy_error_lower_times_1e15"] >= 81600


def test_equivalent_shared_couplings_do_not_trigger_obstruction():
    # Distinct normalizations and a left/right laboratory-axis change.
    j = [[c.iv.mpf(i) for i in row] for row in ((1, 2, 0), (0, -2, 1), (0, 0, 3))]
    k = [[-7 * j[(r + 1) % 3][(s + 2) % 3] for s in range(3)] for r in range(3)]
    assert c.orbit_lower(j, k).a == 0
    with pytest.raises(ValueError, match="shared shape has not been excluded"):
        c.certify_chart([j, k], [j, k], reciprocal=False)


def test_untagged_tomography_and_independent_source_spectrum():
    source = json.loads(c.RADAR.read_text())
    native = json.loads((c.HERE.parents[1] / "results.json").read_text())
    ref = np.asarray(native["x_real"]) + 1j * np.asarray(native["x_imag"])
    a, v = [
        np.subtract(*np.asarray(source[key], dtype=float).mean(axis=-1))
        / source["scale"]
        for key in ("area", "volume")
    ]
    w = a - (np.vdot(a, v) / np.vdot(v, v)) * v
    paulis = [
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    ]
    joint = sum(
        w[i, j] * np.kron(p, q) / 2
        for i, p in enumerate(paulis)
        for j, q in enumerate(paulis)
    )
    d = np.stack((-joint, joint)) / (np.sqrt(2) * np.linalg.norm(w))
    assert np.linalg.norm(d.sum(axis=0)) == 0
    as_two_heads = d.reshape(2, 2, 2, 2, 2)
    for axes in ((1, 3), (2, 4)):
        assert np.linalg.norm(np.trace(as_two_heads, axis1=axes[0], axis2=axes[1])) == 0
    eig, vec = np.linalg.eigh(ref)
    modular = -(vec * np.log(eig)[:, None, :]) @ vec.conj().transpose(0, 2, 1)
    assert 0.006092217 < np.vdot(d, modular).real < 0.006093092
    assert abs(np.vdot(w, v)) < 1e-18
