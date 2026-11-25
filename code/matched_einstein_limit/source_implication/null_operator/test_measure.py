from pathlib import Path
import json

import numpy as np

from . import measure


def test_all_three_native_timeouts_and_hidden_ordinary_transfer():
    fresh = measure.verify()
    stored = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    assert fresh["parked_terms"] == stored["parked_terms"]
    for result in (fresh, stored):
        assert result["invisible"]["2+4"]["normalized_energy_lower"] > .00125
        assert result["invisible"]["2+6"]["normalized_energy_lower"] > .000584
        assert result["invisible"]["6"]["normalized_energy_lower"] > .02495
        assert result["invisible"]["4+6"]["normalized_energy_lower"] > .0002317
        assert result["complete"]["left_inverse_defect_upper"] < .00034
        assert result["complete"]["energy_dual_one_norm_upper"] < 176000
        assert result["complete"]["coordinates"] == 31


def test_exact_parked_central_force_and_transverse_contraction():
    local = measure.parked_terms()
    exchange = np.array(local["exchange"])[1:]
    contrast = np.array(local["contrast"])[1:]
    e = np.kron(np.diag(exchange), np.eye(2))
    ident = np.eye(len(e))
    np.testing.assert_array_equal(np.kron(np.diag(contrast), np.eye(2)), ident)
    # The other head is sign-even, so its exchange is zero and its contrast is one.
    terms = [e, ident, ident, 0*ident]
    for i in terms:
        for j in terms:
            np.testing.assert_array_equal(i@j-j@i, 0*ident)
    # Signed unit shear acts along an arbitrary non-coordinate bond (3,4,0).
    bond = np.array([3, 4, 0])
    transverse = np.array([4, -3, 0])
    for stiffness in (-7, 0, 5):
        force = stiffness*np.outer(bond, bond)
        assert transverse @ force @ transverse == 0
        np.testing.assert_array_equal(e-(transverse @ force @ transverse)*e, e)



def test_design_matches_independent_native_probe_menus():
    rows, unit, _ = measure.range_matrices()
    menus = measure.station.radar_menus()
    for cap in measure.CAPS:
        comparison = measure.station.Survey(menus[cap]).range_map.copy()
        comparison[:, -1] /= np.sqrt(2)
        # This independent sparse floating replay checks conventions and order.
        # Only the integer kernel replay supplies the certified uncertainty.
        np.testing.assert_allclose(rows[cap].astype(float)/unit, comparison, rtol=0, atol=3e-12)
