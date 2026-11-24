import json

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

from . import measure


def test_modular_scaling_with_noncommuting_heads_and_visible_flag():
    rng = np.random.default_rng(7521)
    x = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    joint = x @ x.conj().T
    joint /= np.trace(joint)
    prod = np.kron(*measure.formation.marginals(joint))
    d = joint - prod
    assert np.linalg.norm(d @ prod - prod @ d) > 1e-3
    c, _ = measure.relative_curvature(prod, d)
    def resolvent(t):
        w = np.linalg.solve(prod + t*np.eye(4), d)
        return np.trace(w @ w).real
    np.testing.assert_allclose(c, quad(resolvent, 0, np.inf)[0], rtol=1e-10)
    for row in measure.regional_charge(joint, prod, [])["samples"]:
        assert row["lower"] <= row["scaled_magnitude"]
        if row["upper"] is not None:
            assert row["scaled_magnitude"] <= row["upper"]
    assert measure.charge_slope(prod, d, 0) == 0
    np.testing.assert_allclose(measure.charge_slope(prod, d, 1e-5)/1e-10, -c, rtol=1e-5)
    p = 0.23
    block = np.zeros((8,8), complex)
    block[:4,:4], block[4:,4:] = p*joint, (1-p)*prod
    dx = np.zeros_like(block)
    dx[:4,:4] = p*d
    np.testing.assert_allclose(-np.trace(dx @ measure.log_density(block)),
                               p*measure.charge_slope(prod, d, 1), atol=1e-12)


def test_complete_material_records_and_rate_control():
    bank = json.loads(measure.RECEIPT.read_text())
    for board in bank["boards"]:
        for data in board["regional"]:
            finest = data["range_reads"][-1]
            assert 3.97 < finest["unflagged_doubling_exponent"] < 4.01
            assert 1.99 < finest["flagged_doubling_exponent"] < 2.01
        for probe in board["material"]:
            for pair in probe["preparations"]:
                formed, product = pair["formed"], pair["independent"]
                for key in ("mean_energy", "energy_variance"):
                    np.testing.assert_allclose(formed[key], product[key], atol=2e-12, rtol=0)
                assert formed["energy_variance"] > 1
                for arm in (formed, product):
                    for a in arm["appointments"]:
                        probabilities = np.array([a["relay_miss"], a["return_miss"], a["complete"]])
                        assert probabilities.min() > 0
                        np.testing.assert_allclose(probabilities.sum(), 1, atol=1e-13)
                assert abs(formed["appointments"][-1]["complete"]
                           - product["appointments"][-1]["complete"]) > 1e-5

    # An independent dense evolution retains the head density through the relay.
    from scipy import sparse
    rng = np.random.default_rng(48)
    a = rng.normal(size=(12,12)) + 1j*rng.normal(size=(12,12))
    h = (a + a.conj().T) / 2
    h -= np.trace(h) * np.eye(12) / 12
    source = np.diag([.35,.15,.15,.35])
    record = measure.material_read(sparse.csr_matrix(h), (0, 2), {1: (source, np.eye(4)/4)})
    for item in record[0]["formed"]["appointments"]:
        u = expm(-1j*item["time"]*h)
        incoming = np.zeros((12,12), complex)
        incoming[:4,:4] = source
        density = u @ incoming @ u.conj().T
        relayed = np.zeros_like(density)
        relayed[8:,8:] = density[8:,8:]
        final = u @ relayed @ u.conj().T
        np.testing.assert_allclose(item["complete"], np.trace(final[:4,:4]).real, atol=2e-13)
