import json
from pathlib import Path

import numpy as np

from ..geometry_gain import simplex_area_volume
from . import measure


def test_archived_selected_gradient_hessian_and_finite_comparisons():
    r = measure.run()
    assert [
        (s["bracketed"], s["positive_eta"], s["strict_negative_hessian"])
        for s in r["search"]
    ] == [(532, 447, 0), (200, 200, 60)]
    found = r["selected"]
    np.testing.assert_allclose(found["eta"], 8.540583e-5, rtol=1e-6)
    assert np.prod(found["bracket_determinants"]) < 0
    precise = r["checked"]
    assert precise["tangent_gradient_norm"] < 1e-16
    assert max(precise["eigenvalues"]) < -1.3e-9
    np.testing.assert_allclose(
        precise["eigenvalues"], found["eigenvalues"], rtol=1e-8, atol=1e-17
    )
    for sec in precise["equal_volume_secants"]:
        assert abs(sec["volume_error"]) < 1e-40 and sec["objective_change"] < 0
    c = np.array(r["raw_ranges"])
    x = np.array(found["x"])
    np.testing.assert_allclose(
        simplex_area_volume(c[0, 0] + np.einsum("i,ijk->jk", x, measure.affine(c))),
        precise["area_volume_entropy"][:2],
        atol=2e-14,
    )
    for changed in r["held_at_frozen_gain"]:
        assert abs(changed["ordinary_transfer"] - 1) > 0.1
    stored = json.loads(Path(__file__).with_name("results.json").read_text())
    np.testing.assert_allclose(stored["selected"]["x"], x, atol=2e-11)
